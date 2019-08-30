#include <tgmath.h>
#include "march.h"
#include "table.h"

enum { X, Y, Z };
struct Vec {
    double x;
    double y;
    double z;
};
struct Vec;
static int MarchTetrahedron(struct March *, double *tetr, double *);
static int MarchCubes(struct March *,
		      int (*)(struct March *, double, double, double));
static double offset(double, double);
static void normal(struct March *, double, double, double, double *,
		   double *, double *);

static int
MarchCube1(struct March *q, double x, double y, double z)
{
    double cube[8];
    double a;
    int c, i, e, j, idx, flag;
    double norm[3 * 12];
    double vert[3 * 12];
    double h;
    double *n, *v;

    h = q->spacing;
    for (i = 0; i < 8; i++) {
	cube[i] =
	    q->f(x + Offset[i][0] * h,
		 y + Offset[i][1] * h, z + Offset[i][2] * h, q->fdata);
    }
    idx = 0;
    for (i = 0; i < 8; i++) {
	if (cube[i] <= 0)
	    idx |= 1 << i;
    }
    flag = CubeEdgeFlags[idx];
    if (flag == 0) {
	return 0;
    }
    for (e = 0; e < 12; e++) {
	if (flag & (1 << e)) {
	    n = &norm[3 * e];
	    v = &vert[3 * e];
	    a = offset(cube[Connection[e][0]], cube[Connection[e][1]]);
	    v[X] =
		x + (Offset[Connection[e][0]][0] +
		     a * Direction[e][0]) * h;
	    v[Y] =
		y + (Offset[Connection[e][0]][1] +
		     a * Direction[e][1]) * h;
	    v[Z] =
		z + (Offset[Connection[e][0]][2] +
		     a * Direction[e][2]) * h;
	    normal(q, v[X], v[Y], v[Z], &n[X], &n[Y], &n[Z]);
	}
    }
    for (j = 0; j < 5; j++) {
	if (TriangleConnectionTable[idx][3 * j] < 0)
	    break;
	for (c = 0; c < 3; c++) {
	    i = TriangleConnectionTable[idx][3 * j + c];
	    n = &norm[3 * i];
	    v = &vert[3 * i];
	    q->normal(n[X], n[Y], n[Z], q->cdata);
	    q->vertex(v[X], v[Y], v[Z], q->cdata);
	}
    }
    return 0;
}

static int
MarchCube2(struct March *q, double x, double y, double z)
{
    double cube[8];
    double val[4];
    int i, t, InACube;
    double *p, *te;
    double h;
    double pos[3 * 8];
    double tetr[3 * 4];

    h = q->spacing;
    for (i = 0; i < 8; i++) {
	p = &pos[3 * i];
	p[X] = x + Offset[i][0] * h;
	p[Y] = y + Offset[i][1] * h;
	p[Z] = z + Offset[i][2] * h;
    }
    for (i = 0; i < 8; i++) {
	p = &pos[3 * i];
	cube[i] = q->f(p[X], p[Y], p[Z], q->fdata);
    }
    for (t = 0; t < 6; t++) {
	for (i = 0; i < 4; i++) {
	    InACube = TetrahedronsInACube[t][i];
	    p = &pos[3 * InACube];
	    te = &tetr[3 * i];
	    te[X] = p[X];
	    te[Y] = p[Y];
	    te[Z] = p[Z];
	    val[i] = cube[InACube];
	}
	MarchTetrahedron(q, tetr, val);
    }
    return 0;
}

static int
MarchTetrahedron(struct March *q, double *tetr, double *val)
{
    int e, v0, v1, flag, j, c, i, idx = 0;
    double a, b;
    double *te0, *te1, *v, *n;
    double vert[3 * 6];
    double norm[3 * 6];

    for (i = 0; i < 4; i++) {
	if (val[i] <= 0)
	    idx |= 1 << i;
    }
    flag = TetrahedronEdgeFlags[idx];
    if (flag == 0) {
	return 0;
    }
    for (e = 0; e < 6; e++) {
	if (flag & (1 << e)) {
	    v0 = TetrahedronConnection[e][0];
	    v1 = TetrahedronConnection[e][1];
	    a = offset(val[v0], val[v1]);
	    b = 1.0 - a;
	    v = &vert[3 * e];
	    te0 = &tetr[3 * v0];
	    te1 = &tetr[3 * v1];
	    n = &norm[3 * e];
	    v[X] = b * te0[X] + a * te1[X];
	    v[Y] = b * te0[Y] + a * te1[Y];
	    v[Z] = b * te0[Z] + a * te1[Z];
	    normal(q, v[X], v[Y], v[Z], &n[X], &n[Y], &n[Z]);
	}
    }
    for (j = 0; j < 2; j++) {
	if (TetrahedronTriangles[idx][3 * j] < 0)
	    break;
	for (c = 0; c < 3; c++) {
	    i = TetrahedronTriangles[idx][3 * j + c];
	    v = &vert[3 * i];
	    n = &norm[3 * i];
	    q->normal(n[X], n[Y], n[Z], q->cdata);
	    q->vertex(v[X], v[Y], v[Z], q->cdata);
	}
    }
    return 0;
}

static int
MarchCubes(struct March *q,
	   int (*march0) (struct March *, double, double, double))
{
    enum { X, Y, Z };
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

int
march_tetrahedron(struct March *q)
{
    return MarchCubes(q, MarchCube2);
}

double
offset(double a, double b)
{
    double d;

    d = a - b;
    if (d == 0.0)
	return 0.5;
    return a / d;
}

static double
sq(double x)
{
    return x * x;
}

static void
Normalize(double *u, double *v, double *w)
{
    double len;

    len = sqrt(sq(*u) + sq(*v) + sq(*w));
    if (len != 0.0) {
	*u /= len;
	*v /= len;
	*w /= len;
    }
}

static void
normal(struct March *q, double x, double y, double z, double *u, double *v,
       double *w)
{
    double h;

#define F(x, y, z) (q->f)((x), (y), (z), q->fdata)
    h = q->spacing / 10;
    *u = F(x - h, y, z) - F(x + h, y, z);
    *v = F(x, y - h, z) - F(x, y + h, z);
    *w = F(x, y, z - h) - F(x, y, z + h);
    Normalize(u, v, w);
}
