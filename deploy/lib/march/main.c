#include <tgmath.h>
#include "march.h"
#include "table.h"

enum { X, Y, Z };
static double offset(double, double);

static int
cube(double cube[8], int *pn, double *tri)
{
    double a;
    int c, i, j, idx, flag;
    double vert[3 * 12];
    double *v, *o, *dir;
    int n, k;

    idx = 0;
    for (i = 0; i < 8; i++) {
	if (cube[i] <= 0)
	    idx |= 1 << i;
    }
    if ((flag = CubeEdgeFlags[idx]) == 0) {
	*pn = 0;
	return 0;
    }
    for (i = 0; i < 12; i++) {
	if (flag & (1 << i)) {
	    o = Offset[Connection[i][0]];
	    dir = Direction[i];
	    v = &vert[3 * i];
	    a = offset(cube[Connection[i][0]], cube[Connection[i][1]]);
	    v[X] = o[X] + a * dir[X];
	    v[Y] = o[Y] + a * dir[Y];
	    v[Z] = o[Z] + a * dir[Z];
	}
    }

    n = k = 0;
    for (i = 0; i < 5; i++) {
	if (TriangleConnectionTable[idx][3 * i] < 0)
	    break;
	for (c = 0; c < 3; c++) {
	    j = TriangleConnectionTable[idx][3 * i + c];
	    v = &vert[3 * j];
	    tri[k++] = v[X];
	    tri[k++] = v[Y];
	    tri[k++] = v[Z];
	}
	n++;
    }
    *pn = n;
    return 0;
}

static int
tetrahedron0(double *tetr, double *val, int *pn, double *tri)
{
    int e, v0, v1, flag, j, c, i, idx = 0;
    double a, b;
    double *te0, *te1, *v;
    double vert[3 * 6];
    int n;

    n = *pn;
    for (i = 0; i < 4; i++) {
	if (val[i] <= 0)
	    idx |= 1 << i;
    }
    if ((flag = TetrahedronEdgeFlags[idx]) == 0)
	goto end;
    for (e = 0; e < 6; e++) {
	if (flag & (1 << e)) {
	    v0 = TetrahedronConnection[e][0];
	    v1 = TetrahedronConnection[e][1];
	    a = offset(val[v0], val[v1]);
	    b = 1 - a;
	    v = &vert[3 * e];
	    te0 = &tetr[3 * v0];
	    te1 = &tetr[3 * v1];
	    v[X] = b * te0[X] + a * te1[X];
	    v[Y] = b * te0[Y] + a * te1[Y];
	    v[Z] = b * te0[Z] + a * te1[Z];
	}
    }

    for (j = 0; j < 2; j++) {
	if (TetrahedronTriangles[idx][3 * j] < 0)
	    break;
	for (c = 0; c < 3; c++) {
	    i = TetrahedronTriangles[idx][3 * j + c];
	    v = &vert[3 * i];
	    *tri++ = v[X];
	    *tri++ = v[Y];
	    *tri++ = v[Z];
	}
	n++;
    }
end:
    *pn = n;
    return 0;
}


static int
tetrahedron(double cube[8], int *pn, double *tri)
{
    double val[4];
    int i, t, j;
    double *te, *o;
    double tetr[3 * 4];
    int n;

    n = 0;
    for (t = 0; t < 6; t++) {
	for (i = 0; i < 4; i++) {
	    j = TetrahedronsInACube[t][i];
	    o = Offset[j];
	    te = &tetr[3 * i];
	    te[X] = o[X];
	    te[Y] = o[Y];
	    te[Z] = o[Z];
	    val[i] = cube[j];
	}
	tetrahedron0(tetr, val, &n, &tri[9 * n]);
    }

    *pn = n;
    return 0;
}


static double
offset(double a, double b)
{
    double d;

    d = a - b;
    if (d == 0.0)
	return 0.5;
    return a / d;
}

static void
swap(double *u, int i, int j)
{
    double t;

    t = u[i];
    u[i] = u[j];
    u[j] = t;
}

static int
march(int (*algorithm)(double *, int *, double *), double u[8],
      int *pn, double *tri)
{
    int s;

    swap(u, 2, 3);
    swap(u, 6, 7);
    s = algorithm(u, pn, tri);
    swap(u, 2, 3);
    swap(u, 6, 7);
    return s;
}

int
march_cube(double u[8], int *n, double *tri)
{
    return march(cube, u, n, tri);
}

int
march_tetrahedron(double u[8], int *n, double *tri)
{
    return march(tetrahedron, u, n, tri);
}
