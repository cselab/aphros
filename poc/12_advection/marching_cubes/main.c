#include <stdio.h>
#include <tgmath.h>
#include <GL/glut.h>
#include "table.h"
#include "sample.h"
#include "lib.h"
#include "march.h"
#define	USED(x)		if(x);else{}

struct Vec {
    double x;
    double y;
    double z;
};

/***/
struct March march;
struct March *q = &march;
/***/

struct Sample *sample;
static double
f(double x, double y, double z, void *p)
{
    struct Sample *sample;

    sample = p;
    return sample_f(sample, x, y, z);
}

static int
normal(double x, double y, double z, void *p)
{
    double u, v, w;
    USED(p);
    GetColor(x, y, z, &u, &v, &w);
    glColor3f(u, v, w);
    glNormal3f(x, y, z);
    return 0;
}

static int
vertex(double x, double y, double z, void *p)
{
    USED(p);
    glVertex3f(x, y, z);
    return 0;
}

struct Vec;
static void MarchTetrahedron(struct March*, struct Vec *, double *);
static void MarchingCubes(struct March*, int(*) (struct March*, double, double, double));
static float AmbientGreen[] = { 0.00, 0.25, 0.00, 1.00 };
static float AmbientBlue[] = { 0.00, 0.00, 0.25, 1.00 };
static float DiffuseGreen[] = { 0.00, 0.75, 0.00, 1.00 };
static float DiffuseBlue[] = { 0.00, 0.00, 0.75, 1.00 };
static float SpecularWhite[] = { 1.00, 1.00, 1.00, 1.00 };

int n = 16;
double h;
int PolygonMode = GL_FILL;
int Spin = 0;
int Move = 0;
int Light = 0;

static void
Resize(int Width, int Height)
{
    double Aspect, HalfWorldSize = (1.4142135623730950488016887242097 / 2);

    glViewport(0, 0, Width, Height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    if (Width <= Height) {
	Aspect = (double) Height / (double) Width;
	glOrtho(-HalfWorldSize, HalfWorldSize, -HalfWorldSize * Aspect,
		HalfWorldSize * Aspect, -10 * HalfWorldSize,
		10 * HalfWorldSize);
    } else {
	Aspect = (double) Width / (double) Height;
	glOrtho(-HalfWorldSize * Aspect, HalfWorldSize * Aspect,
		-HalfWorldSize, HalfWorldSize, -10 * HalfWorldSize,
		10 * HalfWorldSize);
    }

    glMatrixMode(GL_MODELVIEW);
}


static int
MarchCube1(struct March *q, double x, double y, double z)
{
    extern int CubeEdgeFlags[256];
    extern int TriangleConnectionTable[256][16];

    double cube[8];
    double Offset;
    int Corner, i, Test, Edge, iTriangle, FlagIndex, EdgeFlags;
    struct Vec EdgeNorm[12];
    struct Vec EdgeVertex[12];


    for (i = 0; i < 8; i++) {
	cube[i] =
	    q->f(x + VertexOffset[i][0] * h,
		 y + VertexOffset[i][1] * h, z + VertexOffset[i][2] * h, q->fdata);
    }


    FlagIndex = 0;
    for (Test = 0; Test < 8; Test++) {
	if (cube[Test] <= 0)
	    FlagIndex |= 1 << Test;
    }


    EdgeFlags = CubeEdgeFlags[FlagIndex];


    if (EdgeFlags == 0) {
	return 0;
    }


    for (Edge = 0; Edge < 12; Edge++) {

	if (EdgeFlags & (1 << Edge)) {
	    Offset = GetOffset(cube[EdgeConnection[Edge][0]],
			       cube[EdgeConnection[Edge][1]]);

	    EdgeVertex[Edge].x =
		x + (VertexOffset[EdgeConnection[Edge][0]][0] +
		     Offset * EdgeDirection[Edge][0]) * h;
	    EdgeVertex[Edge].y =
		y + (VertexOffset[EdgeConnection[Edge][0]][1] +
		     Offset * EdgeDirection[Edge][1]) * h;
	    EdgeVertex[Edge].z =
		z + (VertexOffset[EdgeConnection[Edge][0]][2] +
		     Offset * EdgeDirection[Edge][2]) * h;

	    GetNormal(EdgeVertex[Edge].x, EdgeVertex[Edge].y, EdgeVertex[Edge].z,
		      &EdgeNorm[Edge].x, &EdgeNorm[Edge].y, &EdgeNorm[Edge].z,
		      q->f, q->fdata);

	}
    }



    for (iTriangle = 0; iTriangle < 5; iTriangle++) {
	if (TriangleConnectionTable[FlagIndex][3 * iTriangle] < 0)
	    break;

	for (Corner = 0; Corner < 3; Corner++) {
	    i = TriangleConnectionTable[FlagIndex][3 * iTriangle + Corner];
	    q->normal(EdgeNorm[i].x, EdgeNorm[i].y, EdgeNorm[i].z, q->cdata);
	    q->vertex(EdgeVertex[i].x, EdgeVertex[i].y, EdgeVertex[i].z, q->cdata);
	}
    }
    return 0;
}


static int
MarchCube2(struct March *q, double x, double y, double z)
{
    double cube[8];
    double TetrahedronValue[4];
    int i, Tetrahedron, InACube;
    struct Vec CubePosition[8];
    struct Vec TetrahedronPosition[4];


    for (i = 0; i < 8; i++) {
	CubePosition[i].x = x + VertexOffset[i][0] * (q->spacing);
	CubePosition[i].y = y + VertexOffset[i][1] * (q->spacing);
	CubePosition[i].z = z + VertexOffset[i][2] * (q->spacing);
    }


    for (i = 0; i < 8; i++) {
	cube[i] = q->f(CubePosition[i].x,
		       CubePosition[i].y, CubePosition[i].z, q->fdata);
    }

    for (Tetrahedron = 0; Tetrahedron < 6; Tetrahedron++) {
	for (i = 0; i < 4; i++) {
	    InACube = TetrahedronsInACube[Tetrahedron][i];
	    TetrahedronPosition[i].x = CubePosition[InACube].x;
	    TetrahedronPosition[i].y = CubePosition[InACube].y;
	    TetrahedronPosition[i].z = CubePosition[InACube].z;
	    TetrahedronValue[i] = cube[InACube];
	}
	MarchTetrahedron(q, TetrahedronPosition, TetrahedronValue);
    }
    return 0;
}


static int
(*march0) (struct March*) = march_cube;

static void
Keyboard(unsigned char Key, int i, int j)
{
    USED(i);
    USED(j);
    switch (Key) {
    case 'w':
	{
	    if (PolygonMode == GL_LINE) {
		PolygonMode = GL_FILL;
	    } else {
		PolygonMode = GL_LINE;
	    }
	    glPolygonMode(GL_FRONT_AND_BACK, PolygonMode);
	}
	break;
    case '+':
    case '=':
	{
	    ++n;
	    h = 1.0 / n;
	}
	break;
    case '-':
	{
	    if (n > 1) {
		--n;
		h = 1.0 / n;
	    }
	}
	break;
    case 'c':
	march0 = (march0 == march_cube) ? march_tetrahedron : march_cube;
	break;
    case 's':
	{
	    sample_next(sample);
	}
	break;
    case 'l':
	{
	    if (Light) {
		glDisable(GL_LIGHTING);
	    } else {
		glEnable(GL_LIGHTING);
	    }

	    Light = !Light;
	};
    }
}


static void
Special(int Key, int i, int j)
{
    USED(i);
    USED(j);
    switch (Key) {
    case GLUT_KEY_PAGE_UP:
	sample_inc(sample);
	break;
    case GLUT_KEY_PAGE_DOWN:
	sample_dec(sample);
	break;
    case GLUT_KEY_HOME:
	Spin = !Spin;
	break;
    case GLUT_KEY_END:
	Move = !Move;
	break;
    }
}

static void
Idle(void)
{
    glutPostRedisplay();
}

static void
DrawScene(void)
{
    static double Pitch = 0.0;
    static double Yaw = 0.0;
    static double Time = 0.0;

    march.f = f;
    march.fdata = sample;
    march.size[0] = n;
    march.size[1] = n;
    march.size[2] = n;
    march.spacing = h;
    march.normal = normal;
    march.vertex = vertex;
    march.cdata = NULL;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glPushMatrix();

    if (Spin) {
	Pitch += 4.0;
	Yaw += 2.5;
    }
    if (Move) {
	Time += 0.025;
    }

    sample_time(sample, Time);

    glTranslatef(0.0, 0.0, -1.0);
    glRotatef(-Pitch, 1.0, 0.0, 0.0);
    glRotatef(0.0, 0.0, 1.0, 0.0);
    glRotatef(Yaw, 0.0, 0.0, 1.0);

    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    glColor3f(1.0, 1.0, 1.0);
    glutWireCube(1.0);
    glPopAttrib();


    glPushMatrix();
    glTranslatef(-0.5, -0.5, -0.5);
    glBegin(GL_TRIANGLES);
    march0(q);
    glEnd();
    glPopMatrix();


    glPopMatrix();

    glutSwapBuffers();
}


static void
MarchTetrahedron(struct March *q, struct Vec *TetrahedronPosition, double *TetrahedronValue)
{
    extern int TetrahedronEdgeFlags[16];
    extern int TetrahedronTriangles[16][7];

    int Edge, Vert0, Vert1, EdgeFlags, iTriangle, Corner, i, FlagIndex = 0;
    double Offset, InvOffset;
    struct Vec EdgeVertex[6];
    struct Vec EdgeNorm[6];

    for (i = 0; i < 4; i++) {
	if (TetrahedronValue[i] <= 0)
	    FlagIndex |= 1 << i;
    }


    EdgeFlags = TetrahedronEdgeFlags[FlagIndex];


    if (EdgeFlags == 0) {
	return;
    }


    for (Edge = 0; Edge < 6; Edge++) {

	if (EdgeFlags & (1 << Edge)) {
	    Vert0 = TetrahedronEdgeConnection[Edge][0];
	    Vert1 = TetrahedronEdgeConnection[Edge][1];
	    Offset =
		GetOffset(TetrahedronValue[Vert0],
			  TetrahedronValue[Vert1]);
	    InvOffset = 1.0 - Offset;

	    EdgeVertex[Edge].x =
		InvOffset * TetrahedronPosition[Vert0].x +
		Offset * TetrahedronPosition[Vert1].x;
	    EdgeVertex[Edge].y =
		InvOffset * TetrahedronPosition[Vert0].y +
		Offset * TetrahedronPosition[Vert1].y;
	    EdgeVertex[Edge].z =
		InvOffset * TetrahedronPosition[Vert0].z +
		Offset * TetrahedronPosition[Vert1].z;

	    GetNormal(EdgeVertex[Edge].x, EdgeVertex[Edge].y, EdgeVertex[Edge].z,
		      &EdgeNorm[Edge].x, &EdgeNorm[Edge].y, &EdgeNorm[Edge].z,
		      q->f, q->fdata);
	}
    }

    for (iTriangle = 0; iTriangle < 2; iTriangle++) {
	if (TetrahedronTriangles[FlagIndex][3 * iTriangle] < 0)
	    break;

	for (Corner = 0; Corner < 3; Corner++) {
	    i = TetrahedronTriangles[FlagIndex][3 * iTriangle + Corner];
	    q->normal(EdgeNorm[i].x, EdgeNorm[i].y, EdgeNorm[i].z, q->cdata);
	    q->vertex(EdgeVertex[i].x, EdgeVertex[i].y, EdgeVertex[i].z, q->cdata);
	}
    }
}



static void
MarchingCubes(struct March *q, int (*march0) (struct March*, double, double, double))
{
    enum {X, Y, Z};
    int i, j, k;
    int *size;
    double h;
    size = q->size;
    h = q->spacing;
    for (i = 0; i < size[X]; i++)
	for (j = 0; j < size[Y]; j++)
	    for (k = 0; k < size[Z]; k++)
		march0(q, i * h, j * h, k * h);
}

int
main(int argc, char **argv)
{
    float PropertiesAmbient[] = { 0.50, 0.50, 0.50, 1.00 };
    float PropertiesDiffuse[] = { 0.75, 0.75, 0.75, 1.00 };
    float PropertiesSpecular[] = { 1.00, 1.00, 1.00, 1.00 };

    int Width = 640;
    int Height = 480;

    h = 1.0 / n;
    sample = sample_ini();
    glutInit(&argc, argv);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(Width, Height);
    glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
    glutCreateWindow("Marching Cubes");
    glutDisplayFunc(DrawScene);
    glutIdleFunc(Idle);
    glutReshapeFunc(Resize);
    glutKeyboardFunc(Keyboard);
    glutSpecialFunc(Special);

    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClearDepth(1.0);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK, PolygonMode);

    glLightfv(GL_LIGHT0, GL_AMBIENT, PropertiesAmbient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, PropertiesDiffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, PropertiesSpecular);
    glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0);

    glEnable(GL_LIGHT0);

    glMaterialfv(GL_BACK, GL_AMBIENT, AmbientGreen);
    glMaterialfv(GL_BACK, GL_DIFFUSE, DiffuseGreen);
    glMaterialfv(GL_FRONT, GL_AMBIENT, AmbientBlue);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, DiffuseBlue);
    glMaterialfv(GL_FRONT, GL_SPECULAR, SpecularWhite);
    glMaterialf(GL_FRONT, GL_SHININESS, 25.0);

    Resize(Width, Height);

    PrintHelp();
    glutMainLoop();
}
