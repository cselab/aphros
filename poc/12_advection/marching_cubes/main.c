#include <stdio.h>
#include <tgmath.h>
#include "GL/glut.h"
#include "table.h"
#define	USED(x)		if(x);else{}

struct Vec {
    float x;
    float y;
    float z;
};

//These tables are used so that everything can be done in little loops
// that you can look at all at once rather than in pages and pages of
// unrolled code.

//VertexOffset lists the positions, relative to vertex0, of each of
//the 8 vertices of a cube
static const float VertexOffset[][3] = {
    {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0,
							0.0},
    {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0,
							1.0}
};

//EdgeConnection lists the index of the endpoint vertices for each
//of the 12 edges of the cube
static const int EdgeConnection[][2] = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0},
    {4, 5}, {5, 6}, {6, 7}, {7, 4},
    {0, 4}, {1, 5}, {2, 6}, {3, 7}
};

//EdgeDirection lists the direction vector (vertex1-vertex0) for
//each edge in the cube
static const float EdgeDirection[][3] = {
    {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0}, {0.0, -1.0,
							 0.0},
    {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0}, {0.0, -1.0,
							 0.0},
    {0.0, 0.0, 1.0}, {0.0, 0.0, 1.0}, {0.0, 0.0, 1.0}, {0.0, 0.0,
							1.0}
};

//TetrahedronEdgeConnection lists the index of the endpoint
//vertices for each of the 6 edges of the tetrahedron
static const int TetrahedronEdgeConnection[][2] = {
    {0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}
};

//TetrahedronEdgeConnection lists the index of verticies from a
// cube that made up each of the six tetrahedrons within the cube
static const int TetrahedronsInACube[][4] = {
    {0, 5, 1, 6},
    {0, 1, 2, 6},
    {0, 2, 3, 6},
    {0, 3, 7, 6},
    {0, 7, 4, 6},
    {0, 4, 5, 6},
};

static const float AmbientGreen[] = { 0.00, 0.25, 0.00, 1.00 };
static const float AmbientBlue[] = { 0.00, 0.00, 0.25, 1.00 };
static const float DiffuseGreen[] = { 0.00, 0.75, 0.00, 1.00 };
static const float DiffuseBlue[] = { 0.00, 0.00, 0.75, 1.00 };
static const float SpecularWhite[] = { 1.00, 1.00, 1.00, 1.00 };

int n = 16;
float h;
int PolygonMode = GL_FILL;
float TargetValue = 48.0;
float Time = 0.0;
struct Vec SourcePoint[3];
int Spin = 1;
int Move = 1;
int Light = 1;

void Idle();
void DrawScene();
void Resize(int, int);
void Keyboard(unsigned char Key, int i, int j);
void Special(int Key, int i, int j);

void PrintHelp();
void SetTime(float Time);
float Sample1(float x, float y, float z);
float Sample2(float x, float y, float z);
float Sample3(float x, float y, float z);

float (*Sample) (float x, float y, float z) = Sample1;

void MarchingCubes();
void MarchCube1(float x, float y, float z, float Scale);
void MarchCube2(float x, float y, float z, float Scale);

void (*MarchCube) (float x, float y, float z, float Scale) = MarchCube1;

int
main(int argc, char **argv)
{
    float PropertiesAmbient[] = { 0.50, 0.50, 0.50, 1.00 };
    float PropertiesDiffuse[] = { 0.75, 0.75, 0.75, 1.00 };
    float PropertiesSpecular[] = { 1.00, 1.00, 1.00, 1.00 };

    int Width = 640;
    int Height = 480;

    h = 1.0 / n;

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

void
PrintHelp()
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


void
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

void
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
	    if (Sample == Sample1) {
		Sample = Sample2;
	    } else if (Sample == Sample2) {
		Sample = Sample3;
	    } else {
		Sample = Sample1;
	    }
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


void
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

void
Idle()
{
    glutPostRedisplay();
}

void
DrawScene()
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

    SetTime(Time);

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
    MarchingCubes();
    glEnd();
    glPopMatrix();


    glPopMatrix();

    glutSwapBuffers();
}

//GetOffset finds the approximate point of intersection of the surface
// between two points with the values Value1 and Value2
float
GetOffset(float Value1, float Value2, float ValueDesired)
{
    double Delta = Value2 - Value1;

    if (Delta == 0.0) {
	return 0.5;
    }
    return (ValueDesired - Value1) / Delta;
}


//GetColor generates a color from a given position and normal of a point
void
GetColor(struct Vec *Color, struct Vec *Normal)
{
    float x = Normal->x;
    float y = Normal->y;
    float z = Normal->z;

    Color->x =
	(x > 0.0 ? x : 0.0) + (y < 0.0 ? -0.5 * y : 0.0) + (z <
							    0.0 ? -0.5
							    * z : 0.0);
    Color->y =
	(y > 0.0 ? y : 0.0) + (z < 0.0 ? -0.5 * z : 0.0) + (x <
							    0.0 ? -0.5
							    * x : 0.0);
    Color->z =
	(z > 0.0 ? z : 0.0) + (x < 0.0 ? -0.5 * x : 0.0) + (y <
							    0.0 ? -0.5
							    * y : 0.0);
}

void
NormalizeVector(struct Vec *VectorResult, struct Vec *VectorSource)
{
    float OldLength;
    float Scale;

    OldLength = sqrtf((VectorSource->x * VectorSource->x) +
		      (VectorSource->y * VectorSource->y) +
		      (VectorSource->z * VectorSource->z));

    if (OldLength == 0.0) {
	VectorResult->x = VectorSource->x;
	VectorResult->y = VectorSource->y;
	VectorResult->z = VectorSource->z;
    } else {
	Scale = 1.0 / OldLength;
	VectorResult->x = VectorSource->x * Scale;
	VectorResult->y = VectorSource->y * Scale;
	VectorResult->z = VectorSource->z * Scale;
    }
}


//Generate a sample data set.  Sample1(), Sample2() and Sample3() define three scalar fields whose
// values vary by the x,y and z coordinates and by the Time value set by SetTime()
void
SetTime(float NewTime)
{
    float Offset;
    int SourceNum;

    for (SourceNum = 0; SourceNum < 3; SourceNum++) {
	SourcePoint[SourceNum].x = 0.5;
	SourcePoint[SourceNum].y = 0.5;
	SourcePoint[SourceNum].z = 0.5;
    }

    Time = NewTime;
    Offset = 1.0 + sinf(Time);
    SourcePoint[0].x *= Offset;
    SourcePoint[1].y *= Offset;
    SourcePoint[2].z *= Offset;
}

//Sample1 finds the distance of (x, y, z) from three moving points
float
Sample1(float x, float y, float z)
{
    double Result = 0.0;
    double Dx, Dy, Dz;

    Dx = x - SourcePoint[0].x;
    Dy = y - SourcePoint[0].y;
    Dz = z - SourcePoint[0].z;
    Result += 0.5 / (Dx * Dx + Dy * Dy + Dz * Dz);

    Dx = x - SourcePoint[1].x;
    Dy = y - SourcePoint[1].y;
    Dz = z - SourcePoint[1].z;
    Result += 1.0 / (Dx * Dx + Dy * Dy + Dz * Dz);

    Dx = x - SourcePoint[2].x;
    Dy = y - SourcePoint[2].y;
    Dz = z - SourcePoint[2].z;
    Result += 1.5 / (Dx * Dx + Dy * Dy + Dz * Dz);

    return Result;
}

//Sample2 finds the distance of (x, y, z) from three moving lines
float
Sample2(float x, float y, float z)
{
    double Result = 0.0;
    double Dx, Dy, Dz;

    Dx = x - SourcePoint[0].x;
    Dy = y - SourcePoint[0].y;
    Result += 0.5 / (Dx * Dx + Dy * Dy);

    Dx = x - SourcePoint[1].x;
    Dz = z - SourcePoint[1].z;
    Result += 0.75 / (Dx * Dx + Dz * Dz);

    Dy = y - SourcePoint[2].y;
    Dz = z - SourcePoint[2].z;
    Result += 1.0 / (Dy * Dy + Dz * Dz);

    return Result;
}


//Sample2 defines a height field by plugging the distance from the center into the sin and cos functions
float
Sample3(float x, float y, float z)
{
    float Height =
	20.0 * (Time +
		sqrt((0.5 - x) * (0.5 - x) + (0.5 - y) * (0.5 - y)));
    Height = 1.5 + 0.1 * (sinf(Height) + cosf(Height));
    double Result = (Height - z) * 50.0;

    return Result;
}


//GetNormal() finds the gradient of the scalar field at a point
//This gradient can be used as a very accurate vertx normal for lighting calculations
void
GetNormal(struct Vec *Normal, float x, float y, float z)
{
    Normal->x = Sample(x - 0.01, y, z) - Sample(x + 0.01, y, z);
    Normal->y = Sample(x, y - 0.01, z) - Sample(x, y + 0.01, z);
    Normal->z = Sample(x, y, z - 0.01) - Sample(x, y, z + 0.01);
    NormalizeVector(Normal, Normal);
}


//MarchCube1 performs the Marching Cubes algorithm on a single cube
void
MarchCube1(float x, float y, float z, float Scale)
{
    extern int CubeEdgeFlags[256];
    extern int a2iTriangleConnectionTable[256][16];

    int Corner, i, Test, Edge, iTriangle, FlagIndex, EdgeFlags;
    float Offset;
    struct Vec Color;
    float CubeValue[8];
    struct Vec EdgeVertex[12];
    struct Vec EdgeNorm[12];

    //Make a local copy of the values at the cube's corners
    for (i = 0; i < 8; i++) {
	CubeValue[i] =
	    Sample(x + VertexOffset[i][0] * Scale,
		   y + VertexOffset[i][1] * Scale,
		   z + VertexOffset[i][2] * Scale);
    }

    //Find which vertices are inside of the surface and which are outside
    FlagIndex = 0;
    for (Test = 0; Test < 8; Test++) {
	if (CubeValue[Test] <= TargetValue)
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
	    Offset = GetOffset(CubeValue[EdgeConnection[Edge][0]],
			       CubeValue[EdgeConnection[Edge][1]],
			       TargetValue);

	    EdgeVertex[Edge].x =
		x + (VertexOffset[EdgeConnection[Edge][0]][0] +
		     Offset * EdgeDirection[Edge][0]) * Scale;
	    EdgeVertex[Edge].y =
		y + (VertexOffset[EdgeConnection[Edge][0]][1] +
		     Offset * EdgeDirection[Edge][1]) * Scale;
	    EdgeVertex[Edge].z =
		z + (VertexOffset[EdgeConnection[Edge][0]][2] +
		     Offset * EdgeDirection[Edge][2]) * Scale;

	    GetNormal(&EdgeNorm[Edge], EdgeVertex[Edge].x,
		      EdgeVertex[Edge].y, EdgeVertex[Edge].z);
	}
    }


    //Draw the triangles that were found.  There can be up to five per cube
    for (iTriangle = 0; iTriangle < 5; iTriangle++) {
	if (a2iTriangleConnectionTable[FlagIndex][3 * iTriangle] < 0)
	    break;

	for (Corner = 0; Corner < 3; Corner++) {
	    i = a2iTriangleConnectionTable[FlagIndex][3 * iTriangle +
						      Corner];

	    GetColor(&Color, &EdgeNorm[i]);
	    glColor3f(Color.x, Color.y, Color.z);
	    glNormal3f(EdgeNorm[i].x, EdgeNorm[i].y, EdgeNorm[i].z);
	    glVertex3f(EdgeVertex[i].x, EdgeVertex[i].y, EdgeVertex[i].z);
	}
    }
}

//MarchTetrahedron performs the Marching Tetrahedrons algorithm on a single tetrahedron
void
MarchTetrahedron(struct Vec *TetrahedronPosition, float *TetrahedronValue)
{
    extern int TetrahedronEdgeFlags[16];
    extern int a2iTetrahedronTriangles[16][7];

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
	if (a2iTetrahedronTriangles[FlagIndex][3 * iTriangle] < 0)
	    break;

	for (Corner = 0; Corner < 3; Corner++) {
	    i = a2iTetrahedronTriangles[FlagIndex][3 * iTriangle + Corner];
	    GetColor(&Color, &EdgeNorm[i]);
	    glColor3f(Color.x, Color.y, Color.z);
	    glNormal3f(EdgeNorm[i].x, EdgeNorm[i].y, EdgeNorm[i].z);
	    glVertex3f(EdgeVertex[i].x, EdgeVertex[i].y, EdgeVertex[i].z);
	}
    }
}



//MarchCube2 performs the Marching Tetrahedrons algorithm on a single cube by making six calls to MarchTetrahedron
void
MarchCube2(float x, float y, float z, float Scale)
{
    int i, Tetrahedron, InACube;
    struct Vec CubePosition[8];
    float CubeValue[8];
    struct Vec TetrahedronPosition[4];
    float TetrahedronValue[4];

    //Make a local copy of the cube's corner positions
    for (i = 0; i < 8; i++) {
	CubePosition[i].x = x + VertexOffset[i][0] * Scale;
	CubePosition[i].y = y + VertexOffset[i][1] * Scale;
	CubePosition[i].z = z + VertexOffset[i][2] * Scale;
    }

    //Make a local copy of the cube's corner values
    for (i = 0; i < 8; i++) {
	CubeValue[i] = Sample(CubePosition[i].x,
			      CubePosition[i].y, CubePosition[i].z);
    }

    for (Tetrahedron = 0; Tetrahedron < 6; Tetrahedron++) {
	for (i = 0; i < 4; i++) {
	    InACube = TetrahedronsInACube[Tetrahedron][i];
	    TetrahedronPosition[i].x = CubePosition[InACube].x;
	    TetrahedronPosition[i].y = CubePosition[InACube].y;
	    TetrahedronPosition[i].z = CubePosition[InACube].z;
	    TetrahedronValue[i] = CubeValue[InACube];
	}
	MarchTetrahedron(TetrahedronPosition, TetrahedronValue);
    }
}


//MarchingCubes iterates over the entire dataset, calling MarchCube on each cube
void
MarchingCubes()
{
    int i, j, k;

    for (i = 0; i < n; i++)
	for (j = 0; j < n; j++)
	    for (k = 0; k < n; k++)
		MarchCube(i * h, j * h, k * h, h);
}
