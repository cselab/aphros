#include <stdio.h>
#include "memory.h"
#include "err.h"
#include "bbox.h"

static char me[] = "inside";

enum {
    X, Y, Z
};

struct Bbox {
    double lo[3];
    double hi[3];
};

static double
array_min(int n, const double a[])
{
    int i;
    double x;

    for (i = 0, x = a[0]; i < n; i++)
        if (a[i] < x)
            x = a[i];
    return x;
}

static double
array_max(int n, const double a[])
{
    int i;
    double x;

    for (i = 0, x = a[0]; i < n; i++)
        if (a[i] > x)
            x = a[i];
    return x;
}

#define DI(D, d) \
	do { \
		q->lo[D] = array_min(n, d); \
		q->hi[D] = array_max(n, d); \
	} while (0)

int
bbox_ini(struct Bbox** pq)
{
    struct Bbox*q;

    MALLOC(1, &q);
    *pq = q;
    return 0;
}

int
bbox_fin(struct Bbox* q)
{
    FREE(q);
    return 0;
}

int
bbox_update(struct Bbox* q, int n, const double * x, const double * y, const double * z)
{
    DI(X, x);
    DI(Y, y);
    DI(Z, z);
    return 0;
}

int
bbox_inside(struct Bbox* q, double x, double y, double z)
{
#define CM(d, D) (lo[D] < d && d < hi[D])
    double *lo, *hi;

    lo = q->lo;
    hi = q->hi;
    return CM(x, X) && CM(y, Y) && CM(z, Z);
}

int
bbox_lo(struct Bbox* q, double ** p)
{
    *p = q->lo;
    return 0;
}

int
bbox_hi(struct Bbox* q, double ** p)
{
    *p = q->hi;
    return 0;
}

double
bbox_xhi(struct Bbox* q)
{
    return q->hi[X];
}

double
bbox_zhi(struct Bbox* q)
{
    return q->hi[Z];
}

int
bbox_center(struct Bbox* q, double c[3])
{
    c[X] = (q->lo[X] + q->hi[X]) / 2;
    c[Y] = (q->lo[Y] + q->hi[Y]) / 2;
    c[Z] = (q->lo[Z] + q->hi[Z]) / 2;
    return 0;
}
