#include "march.h"
#include "table.h"
#include "lib.h"

struct Vec {
    double x;
    double y;
    double z;
};
struct Vec;
static int MarchTetrahedron(struct March*, struct Vec *, double *);
static int MarchCubes(struct March*, int(*) (struct March*, double, double, double));

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
    double h;

    h = q->spacing;
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
	return 0;
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
    return 0;
}



static int
MarchCubes(struct March *q, int (*march0) (struct March*, double, double, double))
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
    return 0;
}

int
march_cube(struct March *q)
{
    return MarchCubes(q, MarchCube1);
}
    
int march_tetrahedron(struct March *q)
{
    return MarchCubes(q, MarchCube2);
}

