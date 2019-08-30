#include <tgmath.h>
#include "march.h"
#include "table.h"

enum { X, Y, Z };
static double offset(double, double);

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
MarchCube(struct March *q, double r[3])
{
    double cube[8];
    double a;
    int c, i, j, idx, flag;
    double vert[3 * 12];
    double h;
    double *v, *o, *dir;

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
	    v = &vert[3 * i];
	    a = offset(cube[Connection[i][0]], cube[Connection[i][1]]);
	    v[X] = r[X] + (o[X] + a * dir[X]) * h;
	    v[Y] = r[Y] + (o[Y] + a * dir[Y]) * h;
	    v[Z] = r[Z] + (o[Z] + a * dir[Z]) * h;
	}
    }
    for (i = 0; i < 5; i++) {
	if (TriangleConnectionTable[idx][3 * i] < 0)
	    break;
	for (c = 0; c < 3; c++) {
	    j = TriangleConnectionTable[idx][3 * i + c];
	    v = &vert[3 * j];
	    q->vertex(v, q->cdata);
	}
    }
    return 0;
}

int
march_cube(struct March *q)
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
		MarchCube(q, r);
	    }
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
