#include <stdio.h>
#include <tgmath.h>
#include "GL/glut.h"
#include "table.h"
#define	USED(x)		if(x);else{}

struct Vec {
    float fX;
    float fY;
    float fZ;
};

//These tables are used so that everything can be done in little loops
// that you can look at all at once rather than in pages and pages of
// unrolled code.

//a2fVertexOffset lists the positions, relative to vertex0, of each of
//the 8 vertices of a cube
static const float a2fVertexOffset[][3] = {
    { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 1.0, 1.0, 0.0 }, { 0.0, 1.0,
							      0.0 },
    { 0.0, 0.0, 1.0 }, { 1.0, 0.0, 1.0 }, { 1.0, 1.0, 1.0 }, { 0.0, 1.0,
							      1.0 }
};

//a2iEdgeConnection lists the index of the endpoint vertices for each
//of the 12 edges of the cube
static const int a2iEdgeConnection[][2] = {
    { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 },
    { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 },
    { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 }
};

//a2fEdgeDirection lists the direction vector (vertex1-vertex0) for
//each edge in the cube
static const float a2fEdgeDirection[][3] = {
    { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { -1.0, 0.0, 0.0 }, { 0.0, -1.0,
							       0.0 },
    { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { -1.0, 0.0, 0.0 }, { 0.0, -1.0,
							       0.0 },
    { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0,
							      1.0 }
};

//a2iTetrahedronEdgeConnection lists the index of the endpoint
//vertices for each of the 6 edges of the tetrahedron
static const int a2iTetrahedronEdgeConnection[][2] = {
    { 0, 1 }, { 1, 2 }, { 2, 0 }, { 0, 3 }, { 1, 3 }, { 2, 3 }
};

//a2iTetrahedronEdgeConnection lists the index of verticies from a
// cube that made up each of the six tetrahedrons within the cube
static const int a2iTetrahedronsInACube[][4] = {
    { 0, 5, 1, 6 },
    { 0, 1, 2, 6 },
    { 0, 2, 3, 6 },
    { 0, 3, 7, 6 },
    { 0, 7, 4, 6 },
    { 0, 4, 5, 6 },
};

static const float afAmbientGreen[] = { 0.00, 0.25, 0.00, 1.00 };
static const float afAmbientBlue[] = { 0.00, 0.00, 0.25, 1.00 };
static const float afDiffuseGreen[] = { 0.00, 0.75, 0.00, 1.00 };
static const float afDiffuseBlue[] = { 0.00, 0.00, 0.75, 1.00 };
static const float afSpecularWhite[] = { 1.00, 1.00, 1.00, 1.00 };

int n = 16;
float h;
int ePolygonMode = GL_FILL;
float fTargetValue = 48.0;
float fTime = 0.0;
struct Vec sSourcePoint[3];
int bSpin = 1;
int bMove = 1;
int bLight = 1;

void vIdle();
void vDrawScene();
void vResize(int, int);
void vKeyboard(unsigned char cKey, int i, int j);
void vSpecial(int iKey, int i, int j);

void vPrintHelp();
void vSetTime(float fTime);
float fSample1(float fX, float fY, float fZ);
float fSample2(float fX, float fY, float fZ);
float fSample3(float fX, float fY, float fZ);

float (*fSample)(float fX, float fY, float fZ) = fSample1;

void vMarchingCubes();
void vMarchCube1(float fX, float fY, float fZ, float fScale);
void vMarchCube2(float fX, float fY, float fZ, float fScale);

void (*vMarchCube)(float fX, float fY, float fZ, float fScale) =
    vMarchCube1;

int
main(int argc, char **argv)
{
    float afPropertiesAmbient[] = { 0.50, 0.50, 0.50, 1.00 };
    float afPropertiesDiffuse[] = { 0.75, 0.75, 0.75, 1.00 };
    float afPropertiesSpecular[] = { 1.00, 1.00, 1.00, 1.00 };

    int iWidth = 640.0;
    int iHeight = 480.0;

    h = 1.0 / n;

    glutInit(&argc, argv);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(iWidth, iHeight);
    glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
    glutCreateWindow("Marching Cubes");
    glutDisplayFunc(vDrawScene);
    glutIdleFunc(vIdle);
    glutReshapeFunc(vResize);
    glutKeyboardFunc(vKeyboard);
    glutSpecialFunc(vSpecial);

    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClearDepth(1.0);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK, ePolygonMode);

    glLightfv(GL_LIGHT0, GL_AMBIENT, afPropertiesAmbient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, afPropertiesDiffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, afPropertiesSpecular);
    glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0);

    glEnable(GL_LIGHT0);

    glMaterialfv(GL_BACK, GL_AMBIENT, afAmbientGreen);
    glMaterialfv(GL_BACK, GL_DIFFUSE, afDiffuseGreen);
    glMaterialfv(GL_FRONT, GL_AMBIENT, afAmbientBlue);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, afDiffuseBlue);
    glMaterialfv(GL_FRONT, GL_SPECULAR, afSpecularWhite);
    glMaterialf(GL_FRONT, GL_SHININESS, 25.0);

    vResize(iWidth, iHeight);

    vPrintHelp();
    glutMainLoop();
}

void
vPrintHelp()
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
vResize(int iWidth, int iHeight)
{
    float fAspect, fHalfWorldSize =
	(1.4142135623730950488016887242097 / 2);

    glViewport(0, 0, iWidth, iHeight);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    if (iWidth <= iHeight) {
	fAspect = (float) iHeight / (float) iWidth;
	glOrtho(-fHalfWorldSize, fHalfWorldSize, -fHalfWorldSize * fAspect,
		fHalfWorldSize * fAspect, -10 * fHalfWorldSize,
		10 * fHalfWorldSize);
    } else {
	fAspect = (float) iWidth / (float) iHeight;
	glOrtho(-fHalfWorldSize * fAspect, fHalfWorldSize * fAspect,
		-fHalfWorldSize, fHalfWorldSize, -10 * fHalfWorldSize,
		10 * fHalfWorldSize);
    }

    glMatrixMode(GL_MODELVIEW);
}

void
vKeyboard(unsigned char cKey, int i, int j)
{
    USED(i);
    USED(j);
    switch (cKey) {
    case 'w':
	{
	    if (ePolygonMode == GL_LINE) {
		ePolygonMode = GL_FILL;
	    } else {
		ePolygonMode = GL_LINE;
	    }
	    glPolygonMode(GL_FRONT_AND_BACK, ePolygonMode);
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
	    if (vMarchCube == vMarchCube1) {
		vMarchCube = vMarchCube2;	//Use Marching Tetrahedrons
	    } else {
		vMarchCube = vMarchCube1;	//Use Marching Cubes
	    }
	}
	break;
    case 's':
	{
	    if (fSample == fSample1) {
		fSample = fSample2;
	    } else if (fSample == fSample2) {
		fSample = fSample3;
	    } else {
		fSample = fSample1;
	    }
	}
	break;
    case 'l':
	{
	    if (bLight) {
		glDisable(GL_LIGHTING);	//use vertex colors
	    } else {
		glEnable(GL_LIGHTING);	//use lit material color
	    }

	    bLight = !bLight;
	};
    }
}


void
vSpecial(int iKey, int i, int j)
{
    USED(i);
    USED(j);
    switch (iKey) {
    case GLUT_KEY_PAGE_UP:
	{
	    if (fTargetValue < 1000.0) {
		fTargetValue *= 1.1;
	    }
	}
	break;
    case GLUT_KEY_PAGE_DOWN:
	{
	    if (fTargetValue > 1.0) {
		fTargetValue /= 1.1;
	    }
	}
	break;
    case GLUT_KEY_HOME:
	{
	    bSpin = !bSpin;
	}
	break;
    case GLUT_KEY_END:
	{
	    bMove = !bMove;
	}
	break;
    }
}

void
vIdle()
{
    glutPostRedisplay();
}

void
vDrawScene()
{
    static float fPitch = 0.0;
    static float fYaw = 0.0;
    static float fTime = 0.0;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glPushMatrix();

    if (bSpin) {
	fPitch += 4.0;
	fYaw += 2.5;
    }
    if (bMove) {
	fTime += 0.025;
    }

    vSetTime(fTime);

    glTranslatef(0.0, 0.0, -1.0);
    glRotatef(-fPitch, 1.0, 0.0, 0.0);
    glRotatef(0.0, 0.0, 1.0, 0.0);
    glRotatef(fYaw, 0.0, 0.0, 1.0);

    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    glColor3f(1.0, 1.0, 1.0);
    glutWireCube(1.0);
    glPopAttrib();


    glPushMatrix();
    glTranslatef(-0.5, -0.5, -0.5);
    glBegin(GL_TRIANGLES);
    vMarchingCubes();
    glEnd();
    glPopMatrix();


    glPopMatrix();

    glutSwapBuffers();
}

//fGetOffset finds the approximate point of intersection of the surface
// between two points with the values fValue1 and fValue2
float
fGetOffset(float fValue1, float fValue2, float fValueDesired)
{
    double fDelta = fValue2 - fValue1;

    if (fDelta == 0.0) {
	return 0.5;
    }
    return (fValueDesired - fValue1) / fDelta;
}


//vGetColor generates a color from a given position and normal of a point
void
vGetColor(struct Vec *rfColor, struct Vec *rfNormal)
{
    float fX = rfNormal->fX;
    float fY = rfNormal->fY;
    float fZ = rfNormal->fZ;

    rfColor->fX =
	(fX > 0.0 ? fX : 0.0) + (fY < 0.0 ? -0.5 * fY : 0.0) + (fZ <
								0.0 ? -0.5
								*
								fZ : 0.0);
    rfColor->fY =
	(fY > 0.0 ? fY : 0.0) + (fZ < 0.0 ? -0.5 * fZ : 0.0) + (fX <
								0.0 ? -0.5
								*
								fX : 0.0);
    rfColor->fZ =
	(fZ > 0.0 ? fZ : 0.0) + (fX < 0.0 ? -0.5 * fX : 0.0) + (fY <
								0.0 ? -0.5
								*
								fY : 0.0);
}

void
vNormalizeVector(struct Vec *rfVectorResult, struct Vec *rfVectorSource)
{
    float fOldLength;
    float fScale;

    fOldLength = sqrtf((rfVectorSource->fX * rfVectorSource->fX) +
		       (rfVectorSource->fY * rfVectorSource->fY) +
		       (rfVectorSource->fZ * rfVectorSource->fZ));

    if (fOldLength == 0.0) {
	rfVectorResult->fX = rfVectorSource->fX;
	rfVectorResult->fY = rfVectorSource->fY;
	rfVectorResult->fZ = rfVectorSource->fZ;
    } else {
	fScale = 1.0 / fOldLength;
	rfVectorResult->fX = rfVectorSource->fX * fScale;
	rfVectorResult->fY = rfVectorSource->fY * fScale;
	rfVectorResult->fZ = rfVectorSource->fZ * fScale;
    }
}


//Generate a sample data set.  fSample1(), fSample2() and fSample3() define three scalar fields whose
// values vary by the X,Y and Z coordinates and by the fTime value set by vSetTime()
void
vSetTime(float fNewTime)
{
    float fOffset;
    int iSourceNum;

    for (iSourceNum = 0; iSourceNum < 3; iSourceNum++) {
	sSourcePoint[iSourceNum].fX = 0.5;
	sSourcePoint[iSourceNum].fY = 0.5;
	sSourcePoint[iSourceNum].fZ = 0.5;
    }

    fTime = fNewTime;
    fOffset = 1.0 + sinf(fTime);
    sSourcePoint[0].fX *= fOffset;
    sSourcePoint[1].fY *= fOffset;
    sSourcePoint[2].fZ *= fOffset;
}

//fSample1 finds the distance of (fX, fY, fZ) from three moving points
float
fSample1(float fX, float fY, float fZ)
{
    double fResult = 0.0;
    double fDx, fDy, fDz;

    fDx = fX - sSourcePoint[0].fX;
    fDy = fY - sSourcePoint[0].fY;
    fDz = fZ - sSourcePoint[0].fZ;
    fResult += 0.5 / (fDx * fDx + fDy * fDy + fDz * fDz);

    fDx = fX - sSourcePoint[1].fX;
    fDy = fY - sSourcePoint[1].fY;
    fDz = fZ - sSourcePoint[1].fZ;
    fResult += 1.0 / (fDx * fDx + fDy * fDy + fDz * fDz);

    fDx = fX - sSourcePoint[2].fX;
    fDy = fY - sSourcePoint[2].fY;
    fDz = fZ - sSourcePoint[2].fZ;
    fResult += 1.5 / (fDx * fDx + fDy * fDy + fDz * fDz);

    return fResult;
}

//fSample2 finds the distance of (fX, fY, fZ) from three moving lines
float
fSample2(float fX, float fY, float fZ)
{
    double fResult = 0.0;
    double fDx, fDy, fDz;

    fDx = fX - sSourcePoint[0].fX;
    fDy = fY - sSourcePoint[0].fY;
    fResult += 0.5 / (fDx * fDx + fDy * fDy);

    fDx = fX - sSourcePoint[1].fX;
    fDz = fZ - sSourcePoint[1].fZ;
    fResult += 0.75 / (fDx * fDx + fDz * fDz);

    fDy = fY - sSourcePoint[2].fY;
    fDz = fZ - sSourcePoint[2].fZ;
    fResult += 1.0 / (fDy * fDy + fDz * fDz);

    return fResult;
}


//fSample2 defines a height field by plugging the distance from the center into the sin and cos functions
float
fSample3(float fX, float fY, float fZ)
{
    float fHeight =
	20.0 * (fTime +
		sqrt((0.5 - fX) * (0.5 - fX) + (0.5 - fY) * (0.5 - fY)));
    fHeight = 1.5 + 0.1 * (sinf(fHeight) + cosf(fHeight));
    double fResult = (fHeight - fZ) * 50.0;

    return fResult;
}


//vGetNormal() finds the gradient of the scalar field at a point
//This gradient can be used as a very accurate vertx normal for lighting calculations
void
vGetNormal(struct Vec *rfNormal, float fX, float fY, float fZ)
{
    rfNormal->fX = fSample(fX - 0.01, fY, fZ) - fSample(fX + 0.01, fY, fZ);
    rfNormal->fY = fSample(fX, fY - 0.01, fZ) - fSample(fX, fY + 0.01, fZ);
    rfNormal->fZ = fSample(fX, fY, fZ - 0.01) - fSample(fX, fY, fZ + 0.01);
    vNormalizeVector(rfNormal, rfNormal);
}


//vMarchCube1 performs the Marching Cubes algorithm on a single cube
void
vMarchCube1(float fX, float fY, float fZ, float fScale)
{
    extern int aiCubeEdgeFlags[256];
    extern int a2iTriangleConnectionTable[256][16];

    int iCorner, i, iTest, iEdge, iTriangle, iFlagIndex, iEdgeFlags;
    float fOffset;
    struct Vec sColor;
    float afCubeValue[8];
    struct Vec asEdgeVertex[12];
    struct Vec asEdgeNorm[12];

    //Make a local copy of the values at the cube's corners
    for (i = 0; i < 8; i++) {
	afCubeValue[i] =
	    fSample(fX + a2fVertexOffset[i][0] * fScale,
		    fY + a2fVertexOffset[i][1] * fScale,
		    fZ + a2fVertexOffset[i][2] * fScale);
    }

    //Find which vertices are inside of the surface and which are outside
    iFlagIndex = 0;
    for (iTest = 0; iTest < 8; iTest++) {
	if (afCubeValue[iTest] <= fTargetValue)
	    iFlagIndex |= 1 << iTest;
    }

    //Find which edges are intersected by the surface
    iEdgeFlags = aiCubeEdgeFlags[iFlagIndex];

    //If the cube is entirely inside or outside of the surface, then there will be no intersections
    if (iEdgeFlags == 0) {
	return;
    }
    //Find the point of intersection of the surface with each edge
    //Then find the normal to the surface at those points
    for (iEdge = 0; iEdge < 12; iEdge++) {
	//if there is an intersection on this edge
	if (iEdgeFlags & (1 << iEdge)) {
	    fOffset = fGetOffset(afCubeValue[a2iEdgeConnection[iEdge][0]],
				 afCubeValue[a2iEdgeConnection[iEdge][1]],
				 fTargetValue);

	    asEdgeVertex[iEdge].fX =
		fX + (a2fVertexOffset[a2iEdgeConnection[iEdge][0]][0] +
		      fOffset * a2fEdgeDirection[iEdge][0]) * fScale;
	    asEdgeVertex[iEdge].fY =
		fY + (a2fVertexOffset[a2iEdgeConnection[iEdge][0]][1] +
		      fOffset * a2fEdgeDirection[iEdge][1]) * fScale;
	    asEdgeVertex[iEdge].fZ =
		fZ + (a2fVertexOffset[a2iEdgeConnection[iEdge][0]][2] +
		      fOffset * a2fEdgeDirection[iEdge][2]) * fScale;

	    vGetNormal(&asEdgeNorm[iEdge], asEdgeVertex[iEdge].fX,
		       asEdgeVertex[iEdge].fY, asEdgeVertex[iEdge].fZ);
	}
    }


    //Draw the triangles that were found.  There can be up to five per cube
    for (iTriangle = 0; iTriangle < 5; iTriangle++) {
	if (a2iTriangleConnectionTable[iFlagIndex][3 * iTriangle] < 0)
	    break;

	for (iCorner = 0; iCorner < 3; iCorner++) {
	    i = a2iTriangleConnectionTable[iFlagIndex][3 * iTriangle +
						       iCorner];

	    vGetColor(&sColor, &asEdgeNorm[i]);
	    glColor3f(sColor.fX, sColor.fY, sColor.fZ);
	    glNormal3f(asEdgeNorm[i].fX, asEdgeNorm[i].fY,
		       asEdgeNorm[i].fZ);
	    glVertex3f(asEdgeVertex[i].fX, asEdgeVertex[i].fY,
		       asEdgeVertex[i].fZ);
	}
    }
}

//vMarchTetrahedron performs the Marching Tetrahedrons algorithm on a single tetrahedron
void
vMarchTetrahedron(struct Vec *pasTetrahedronPosition,
		  float *pafTetrahedronValue)
{
    extern int aiTetrahedronEdgeFlags[16];
    extern int a2iTetrahedronTriangles[16][7];

    int iEdge, iVert0, iVert1, iEdgeFlags, iTriangle, iCorner, i,
	iFlagIndex = 0;
    float fOffset, fInvOffset;
    struct Vec asEdgeVertex[6];
    struct Vec asEdgeNorm[6];
    struct Vec sColor;

    //Find which vertices are inside of the surface and which are outside
    for (i = 0; i < 4; i++) {
	if (pafTetrahedronValue[i] <= fTargetValue)
	    iFlagIndex |= 1 << i;
    }

    //Find which edges are intersected by the surface
    iEdgeFlags = aiTetrahedronEdgeFlags[iFlagIndex];

    //If the tetrahedron is entirely inside or outside of the surface, then there will be no intersections
    if (iEdgeFlags == 0) {
	return;
    }
    //Find the point of intersection of the surface with each edge
    // Then find the normal to the surface at those points
    for (iEdge = 0; iEdge < 6; iEdge++) {
	//if there is an intersection on this edge
	if (iEdgeFlags & (1 << iEdge)) {
	    iVert0 = a2iTetrahedronEdgeConnection[iEdge][0];
	    iVert1 = a2iTetrahedronEdgeConnection[iEdge][1];
	    fOffset =
		fGetOffset(pafTetrahedronValue[iVert0],
			   pafTetrahedronValue[iVert1], fTargetValue);
	    fInvOffset = 1.0 - fOffset;

	    asEdgeVertex[iEdge].fX =
		fInvOffset * pasTetrahedronPosition[iVert0].fX +
		fOffset * pasTetrahedronPosition[iVert1].fX;
	    asEdgeVertex[iEdge].fY =
		fInvOffset * pasTetrahedronPosition[iVert0].fY +
		fOffset * pasTetrahedronPosition[iVert1].fY;
	    asEdgeVertex[iEdge].fZ =
		fInvOffset * pasTetrahedronPosition[iVert0].fZ +
		fOffset * pasTetrahedronPosition[iVert1].fZ;

	    vGetNormal(&asEdgeNorm[iEdge], asEdgeVertex[iEdge].fX,
		       asEdgeVertex[iEdge].fY, asEdgeVertex[iEdge].fZ);
	}
    }
    //Draw the triangles that were found.  There can be up to 2 per tetrahedron
    for (iTriangle = 0; iTriangle < 2; iTriangle++) {
	if (a2iTetrahedronTriangles[iFlagIndex][3 * iTriangle] < 0)
	    break;

	for (iCorner = 0; iCorner < 3; iCorner++) {
	    i = a2iTetrahedronTriangles[iFlagIndex][3 * iTriangle +
						    iCorner];
	    vGetColor(&sColor, &asEdgeNorm[i]);
	    glColor3f(sColor.fX, sColor.fY, sColor.fZ);
	    glNormal3f(asEdgeNorm[i].fX, asEdgeNorm[i].fY,
		       asEdgeNorm[i].fZ);
	    glVertex3f(asEdgeVertex[i].fX, asEdgeVertex[i].fY,
		       asEdgeVertex[i].fZ);
	}
    }
}



//vMarchCube2 performs the Marching Tetrahedrons algorithm on a single cube by making six calls to vMarchTetrahedron
void
vMarchCube2(float fX, float fY, float fZ, float fScale)
{
    int i, iTetrahedron, iInACube;
    struct Vec asCubePosition[8];
    float afCubeValue[8];
    struct Vec asTetrahedronPosition[4];
    float afTetrahedronValue[4];

    //Make a local copy of the cube's corner positions
    for (i = 0; i < 8; i++) {
	asCubePosition[i].fX = fX + a2fVertexOffset[i][0] * fScale;
	asCubePosition[i].fY = fY + a2fVertexOffset[i][1] * fScale;
	asCubePosition[i].fZ = fZ + a2fVertexOffset[i][2] * fScale;
    }

    //Make a local copy of the cube's corner values
    for (i = 0; i < 8; i++) {
	afCubeValue[i] = fSample(asCubePosition[i].fX,
				 asCubePosition[i].fY,
				 asCubePosition[i].fZ);
    }

    for (iTetrahedron = 0; iTetrahedron < 6; iTetrahedron++) {
	for (i = 0; i < 4; i++) {
	    iInACube = a2iTetrahedronsInACube[iTetrahedron][i];
	    asTetrahedronPosition[i].fX = asCubePosition[iInACube].fX;
	    asTetrahedronPosition[i].fY = asCubePosition[iInACube].fY;
	    asTetrahedronPosition[i].fZ = asCubePosition[iInACube].fZ;
	    afTetrahedronValue[i] = afCubeValue[iInACube];
	}
	vMarchTetrahedron(asTetrahedronPosition, afTetrahedronValue);
    }
}


//vMarchingCubes iterates over the entire dataset, calling vMarchCube on each cube
void
vMarchingCubes()
{
    int i, j, k;

    for (i = 0; i < n; i++)
	for (j = 0; j < n; j++)
	    for (k = 0; k < n; k++)
		vMarchCube(i * h, j * h, k * h, h);
}
