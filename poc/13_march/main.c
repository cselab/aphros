#include <tgmath.h>
#include "march.h"
#include "table.h"

enum { X, Y, Z };
static int MarchTetrahedron(struct March *, double *tetr, double *);
static int MarchCubes(struct March *,
		      int (*)(struct March *, double[3]));
static double offset(double, double);
static void normal(struct March *, double[3], double[3]);

static void
get_cube(struct March *q, double r[3], double cube[8])
{
    int i;
    double *o, p[3];
    double h;

    h = q->spacing;
    for (i = 0; i < 8; i++) {
	o = Offset[i];
	p[X] = r[X] + o[X] * h;
	p[Y] = r[Y] + o[Y] * h;
	p[Z] = r[Z] + o[Z] * h;
	cube[i] =
	    q->f(p, q->fdata);
    }
}

static int
MarchCube1(struct March *q, double r[3])
{
    double cube[8];
    double a;
    int c, i, j, idx, flag;
    double norm[3 * 12];
    double vert[3 * 12];
    double h;
    double *n, *v, *o, *dir;

    h = q->spacing;
    get_cube(q, r, cube);
    
    idx = 0;
    for (i = 0; i < 8; i++) {
	if (cube[i] <= 0)
	    idx |= 1 << i;
    }
    if ((flag = CubeEdgeFlags[idx]) == 0)
	return 0;
    for (i = 0; i < 12; i++) {
	if (flag & (1 << i)) {
	    o = Offset[Connection[i][0]];
	    dir = Direction[i];
	    n = &norm[3 * i];
	    v = &vert[3 * i];
	    a = offset(cube[Connection[i][0]], cube[Connection[i][1]]);
	    v[X] = r[X] + (o[X] + a * dir[X]) * h;
	    v[Y] = r[Y] + (o[Y] + a * dir[Y]) * h;
	    v[Z] = r[Z] + (o[Z] + a * dir[Z]) * h;
	    normal(q, v, n);
	}
    }
    for (i = 0; i < 5; i++) {
	if (TriangleConnectionTable[idx][3 * i] < 0)
	    break;
	for (c = 0; c < 3; c++) {
	    j = TriangleConnectionTable[idx][3 * i + c];
	    n = &norm[3 * j];
	    v = &vert[3 * j];
	    q->normal(n, q->cdata);
	    q->vertex(v, q->cdata);
	}
    }
    return 0;
}

static int
MarchCube2(struct March *q, double r[3])
{
    double cube[8];
    double val[4];
    int i, t, j;
    double *p, *te, *o;
    double h;
    double pos[3 * 8];
    double tetr[3 * 4];

    h = q->spacing;
    for (i = 0; i < 8; i++) {
	p = &pos[3 * i];
	o = Offset[i];
	p[X] = r[X] + o[X] * h;
	p[Y] = r[Y] + o[Y] * h;
	p[Z] = r[Z] + o[Z] * h;
    }
    get_cube(q, r, cube);
    for (t = 0; t < 6; t++) {
	for (i = 0; i < 4; i++) {
	    j = TetrahedronsInACube[t][i];
	    p = &pos[3 * j];
	    te = &tetr[3 * i];
	    te[X] = p[X];
	    te[Y] = p[Y];
	    te[Z] = p[Z];
	    val[i] = cube[j];
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
	    normal(q, v, n);
	}
    }
    for (j = 0; j < 2; j++) {
	if (TetrahedronTriangles[idx][3 * j] < 0)
	    break;
	for (c = 0; c < 3; c++) {
	    i = TetrahedronTriangles[idx][3 * j + c];
	    v = &vert[3 * i];
	    n = &norm[3 * i];
	    q->normal(n, q->cdata);
	    q->vertex(v, q->cdata);
	}
    }
    return 0;
}

static int
MarchCubes(struct March *q,
	   int (*march0) (struct March *, double[3]))
{
    enum { X, Y, Z };
    int i, j, k;
    int *size;
    double h, r[3];

    size = q->size;
    h = q->spacing;
    for (i = 0; i < size[X]; i++)
	for (j = 0; j < size[Y]; j++)
	    for (k = 0; k < size[Z]; k++) {
		r[X] = i * h;
		r[Y] = j * h;
		r[Z] = k * h;
		march0(q, r);
	    }
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

static double
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
Normalize(double v[3])
{
    double len;

    len = sqrt(sq(v[X]) + sq(v[Y]) + sq(v[Z]));
    if (len != 0.0) {
	v[X] /= len;
	v[Y] /= len;
	v[Z] /= len;
    }
}

static double
F(struct March *q, double x, double y, double z)
{
    double p[3];
    p[X] = x;
    p[Y] = y;
    p[Z] = z;
    return (q->f)(p, q->fdata);
}

static void
normal(struct March *q, double r[3], double n[3])
{
    double h, x, y, z;
    h = q->spacing / 10;
    x = r[X];
    y = r[Y];
    z = r[Z];
    n[X] = F(q, x - h, y, z) - F(q, x + h, y, z);
    n[Y] = F(q, x, y - h, z) - F(q, x, y + h, z);
    n[Z] = F(q, x, y, z - h) - F(q, x, y, z + h);
    Normalize(n);
}
