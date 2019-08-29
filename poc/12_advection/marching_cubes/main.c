#include <stdio.h>
#include <tgmath.h>
#include "GL/glut.h"
#include "table.h"
#include "sample.h"
#include "lib.h"
#define	USED(x)		if(x);else{}

struct Sample *sample;
static double f(double x, double y, double z, void *p) {
    struct Sample *sample;
    sample = p;
    return sample_f(sample, x, y, z);
}

struct Vec;
static void GetNormal(struct Vec*, float, float, float);
static void MarchTetrahedron(struct Vec *, float *);
static void MarchingCubes(double(*)(double, double, double, void*), void *);
static float AmbientGreen[] = { 0.00, 0.25, 0.00, 1.00 };
static float AmbientBlue[] = { 0.00, 0.00, 0.25, 1.00 };
static float DiffuseGreen[] = { 0.00, 0.75, 0.00, 1.00 };
static float DiffuseBlue[] = { 0.00, 0.00, 0.75, 1.00 };
static float SpecularWhite[] = { 1.00, 1.00, 1.00, 1.00 };

int n = 16;
float h;
int PolygonMode = GL_FILL;
float TargetValue = 48.0;
int Spin = 1;
int Move = 1;
int Light = 1;

static void
PrintHelp(void)
{
    printf
	("Marching Cubes Example by Cory Bloyd (dejaspaminacan@my-deja.com)\n\n");

    printf("+/-  increase/decrease sample density\n");
    printf("PageUp/PageDown  increase/decrease surface value\n");
    printf("s  change sample function\n");
    printf("c  toggle marching cubes / marching tetrahedrons\n");
    printf("w  wireframe on/off\n");
    printf("l  toggle lighting / color-by-normal\n");
    printf("Home  spin scene on/off\n");
    printf("End  source point animation on/off\n");
}


static void
Resize(int Width, int Height)
{
    float Aspect, HalfWorldSize = (1.4142135623730950488016887242097 / 2);

    glViewport(0, 0, Width, Height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    if (Width <= Height) {
	Aspect = (float) Height / (float) Width;
	glOrtho(-HalfWorldSize, HalfWorldSize, -HalfWorldSize * Aspect,
		HalfWorldSize * Aspect, -10 * HalfWorldSize,
		10 * HalfWorldSize);
    } else {
	Aspect = (float) Width / (float) Height;
	glOrtho(-HalfWorldSize * Aspect, HalfWorldSize * Aspect,
		-HalfWorldSize, HalfWorldSize, -10 * HalfWorldSize,
		10 * HalfWorldSize);
    }

    glMatrixMode(GL_MODELVIEW);
}

//MarchCube1 performs the Marching Cubes algorithm on a single cube
static void
MarchCube1(float x, float y, float z, float h,
	   double f(double, double, double, void*),
	   void *p)
{
    extern int CubeEdgeFlags[256];
    extern int TriangleConnectionTable[256][16];

    float cube[8];
    float Offset;
    int Corner, i, Test, Edge, iTriangle, FlagIndex, EdgeFlags;
    struct Vec Color;
    struct Vec EdgeNorm[12];
    struct Vec EdgeVertex[12];

    //Make a local copy of the values at the cube's corners
    for (i = 0; i < 8; i++) {
	cube[i] =
	    f(x + VertexOffset[i][0] * h,
	      y + VertexOffset[i][1] * h,
	      z + VertexOffset[i][2] * h,
	      p);
    }

    //Find which vertices are inside of the surface and which are outside
    FlagIndex = 0;
    for (Test = 0; Test < 8; Test++) {
	if (cube[Test] <= TargetValue)
	    FlagIndex |= 1 << Test;
    }

    //Find which edges are intersected by the surface
    EdgeFlags = CubeEdgeFlags[FlagIndex];

    //If the cube is entirely inside or outside of the surface, then there will be no intersections
    if (EdgeFlags == 0) {
	return;
    }
    //Find the point of intersection of the surface with each edge
    //Then find the normal to the surface at those points
    for (Edge = 0; Edge < 12; Edge++) {
	//if there is an intersection on this edge
	if (EdgeFlags & (1 << Edge)) {
	    Offset = GetOffset(cube[EdgeConnection[Edge][0]],
			       cube[EdgeConnection[Edge][1]],
			       TargetValue);

	    EdgeVertex[Edge].x =
		x + (VertexOffset[EdgeConnection[Edge][0]][0] +
		     Offset * EdgeDirection[Edge][0]) * h;
	    EdgeVertex[Edge].y =
		y + (VertexOffset[EdgeConnection[Edge][0]][1] +
		     Offset * EdgeDirection[Edge][1]) * h;
	    EdgeVertex[Edge].z =
		z + (VertexOffset[EdgeConnection[Edge][0]][2] +
		     Offset * EdgeDirection[Edge][2]) * h;

	    GetNormal(&EdgeNorm[Edge], EdgeVertex[Edge].x,
		      EdgeVertex[Edge].y, EdgeVertex[Edge].z);
	}
    }


    //Draw the triangles that were found.  There can be up to five per cube
    for (iTriangle = 0; iTriangle < 5; iTriangle++) {
	if (TriangleConnectionTable[FlagIndex][3 * iTriangle] < 0)
	    break;

	for (Corner = 0; Corner < 3; Corner++) {
	    i = TriangleConnectionTable[FlagIndex][3 * iTriangle + Corner];

	    GetColor(&Color, &EdgeNorm[i]);
	    glColor3f(Color.x, Color.y, Color.z);
	    glNormal3f(EdgeNorm[i].x, EdgeNorm[i].y, EdgeNorm[i].z);
	    glVertex3f(EdgeVertex[i].x, EdgeVertex[i].y, EdgeVertex[i].z);
	}
    }
}

//MarchCube2 performs the Marching Tetrahedrons algorithm on a single cube by making six calls to MarchTetrahedron
static void
MarchCube2(float x, float y, float z, float Scale,
	   double f(double, double, double, void*),
	   void *p)
{
    float cube[8];
    float TetrahedronValue[4];
    int i, Tetrahedron, InACube;
    struct Vec CubePosition[8];
    struct Vec TetrahedronPosition[4];

    //Make a local copy of the cube's corner positions
    for (i = 0; i < 8; i++) {
	CubePosition[i].x = x + VertexOffset[i][0] * Scale;
	CubePosition[i].y = y + VertexOffset[i][1] * Scale;
	CubePosition[i].z = z + VertexOffset[i][2] * Scale;
    }

    //Make a local copy of the cube's corner values
    for (i = 0; i < 8; i++) {
	cube[i] = f(CubePosition[i].x,
		    CubePosition[i].y, CubePosition[i].z, p);
    }

    for (Tetrahedron = 0; Tetrahedron < 6; Tetrahedron++) {
	for (i = 0; i < 4; i++) {
	    InACube = TetrahedronsInACube[Tetrahedron][i];
	    TetrahedronPosition[i].x = CubePosition[InACube].x;
	    TetrahedronPosition[i].y = CubePosition[InACube].y;
	    TetrahedronPosition[i].z = CubePosition[InACube].z;
	    TetrahedronValue[i] = cube[InACube];
	}
	MarchTetrahedron(TetrahedronPosition, TetrahedronValue);
    }
}

static void (*MarchCube) (float, float, float, float h,
			  double (*)(double, double, double, void*), void*) = MarchCube1;
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
	{
	    if (MarchCube == MarchCube1) {
		MarchCube = MarchCube2;	//Use Marching Tetrahedrons
	    } else {
		MarchCube = MarchCube1;	//Use Marching Cubes
	    }
	}
	break;
    case 's':
	{
	    sample_next(sample);
	}
	break;
    case 'l':
	{
	    if (Light) {
		glDisable(GL_LIGHTING);	//use vertex colors
	    } else {
		glEnable(GL_LIGHTING);	//use lit material color
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
	{
	    if (TargetValue < 1000.0) {
		TargetValue *= 1.1;
	    }
	}
	break;
    case GLUT_KEY_PAGE_DOWN:
	{
	    if (TargetValue > 1.0) {
		TargetValue /= 1.1;
	    }
	}
	break;
    case GLUT_KEY_HOME:
	{
	    Spin = !Spin;
	}
	break;
    case GLUT_KEY_END:
	{
	    Move = !Move;
	}
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
    static float Pitch = 0.0;
    static float Yaw = 0.0;
    static float Time = 0.0;

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
    MarchingCubes(f, sample);
    glEnd();
    glPopMatrix();


    glPopMatrix();

    glutSwapBuffers();
}

static void
Normalize(struct Vec *v)
{
    float len;

    len = sqrt((v->x * v->x) + (v->y * v->y) + (v->z * v->z));
    if (len != 0.0) {
	v->x /= len;
	v->y /= len;
	v->z /= len;
    }
}

//GetNormal() finds the gradient of the scalar field at a point
//This gradient can be used as a very accurate vertx normal for lighting calculations
static void
GetNormal(struct Vec *Normal, float x, float y, float z)
{
#define F(x, y, z) sample_f(sample, (x), (y), (z))
    Normal->x = F(x - 0.01, y, z) - F(x + 0.01, y, z);
    Normal->y = F(x, y - 0.01, z) - F(x, y + 0.01, z);
    Normal->z = F(x, y, z - 0.01) - F(x, y, z + 0.01);
    Normalize(Normal);
#undef F
}


//MarchTetrahedron performs the Marching Tetrahedrons algorithm on a single tetrahedron
static void
MarchTetrahedron(struct Vec *TetrahedronPosition, float *TetrahedronValue)
{
    extern int TetrahedronEdgeFlags[16];
    extern int TetrahedronTriangles[16][7];

    int Edge, Vert0, Vert1, EdgeFlags, iTriangle, Corner, i, FlagIndex = 0;
    float Offset, InvOffset;
    struct Vec EdgeVertex[6];
    struct Vec EdgeNorm[6];
    struct Vec Color;

    //Find which vertices are inside of the surface and which are outside
    for (i = 0; i < 4; i++) {
	if (TetrahedronValue[i] <= TargetValue)
	    FlagIndex |= 1 << i;
    }

    //Find which edges are intersected by the surface
    EdgeFlags = TetrahedronEdgeFlags[FlagIndex];

    //If the tetrahedron is entirely inside or outside of the surface, then there will be no intersections
    if (EdgeFlags == 0) {
	return;
    }
    //Find the point of intersection of the surface with each edge
    // Then find the normal to the surface at those points
    for (Edge = 0; Edge < 6; Edge++) {
	//if there is an intersection on this edge
	if (EdgeFlags & (1 << Edge)) {
	    Vert0 = TetrahedronEdgeConnection[Edge][0];
	    Vert1 = TetrahedronEdgeConnection[Edge][1];
	    Offset =
		GetOffset(TetrahedronValue[Vert0],
			  TetrahedronValue[Vert1], TargetValue);
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

	    GetNormal(&EdgeNorm[Edge], EdgeVertex[Edge].x,
		      EdgeVertex[Edge].y, EdgeVertex[Edge].z);
	}
    }
    //Draw the triangles that were found.  There can be up to 2 per tetrahedron
    for (iTriangle = 0; iTriangle < 2; iTriangle++) {
	if (TetrahedronTriangles[FlagIndex][3 * iTriangle] < 0)
	    break;

	for (Corner = 0; Corner < 3; Corner++) {
	    i = TetrahedronTriangles[FlagIndex][3 * iTriangle + Corner];
	    GetColor(&Color, &EdgeNorm[i]);
	    glColor3f(Color.x, Color.y, Color.z);
	    glNormal3f(EdgeNorm[i].x, EdgeNorm[i].y, EdgeNorm[i].z);
	    glVertex3f(EdgeVertex[i].x, EdgeVertex[i].y, EdgeVertex[i].z);
	}
    }
}


//MarchingCubes iterates over the entire dataset, calling MarchCube on each cube
static void
MarchingCubes(double f(double, double, double, void*), void *p)
{
    int i, j, k;

    for (i = 0; i < n; i++)
	for (j = 0; j < n; j++)
	    for (k = 0; k < n; k++)
		MarchCube(i * h, j * h, k * h, h, f, p);
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
