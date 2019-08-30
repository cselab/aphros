#include <tgmath.h>
#include "march.h"
#include "table.h"

enum { X, Y, Z };
static double offset(double, double);

int
march_cube(double cube[8], int *pn, double *tri)
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
    if ((flag = CubeEdgeFlags[idx]) == 0)
	return 0;
    for (i = 0; i < 12; i++) {
	if (flag & (1 << i)) {
	    o = Offset[Connection[i][0]];
	    dir = Direction[i];
	    v = &vert[3 * i];
	    a = offset(cube[Connection[i][0]], cube[Connection[i][1]]);
	    v[X] = (o[X] + a * dir[X]);
	    v[Y] = (o[Y] + a * dir[Y]);
	    v[Z] = (o[Z] + a * dir[Z]);
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
	    n++;
	}
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
