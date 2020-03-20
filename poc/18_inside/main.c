#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "err.h"
#include "memory.h"
#include "predicate.h"
#include "inside.h"

static char me[] = "inside";


enum {
    X, Y, Z
};

struct Inside {
    struct He *he;
    const double *x, *y, *z;
    int Update;
};

int
surface_ini(double lo[2], double hi[2], double size, struct Inside ** pq)
{
    struct Inside *q;

    MALLOC(1, &q);
    q->Update = 0;
    //predicate_ini();
    *pq = q;
    return 0;
}

int
surface_fin(struct Inside * q)
{
    FREE(q);
    return 0;
}

int
surface_update(struct Inside * q, struct He * he, const double * x, const double * y,
               const double * z)
{
    int n;

    q->he = he;
    q->x = x;
    q->y = y;
    q->z = z;
    q->Update = 1;
    return 0;
}

#define max(a, b) ( (a) > (b) ? (a) : (b) )
int
surface_inside(struct Inside * q, double u, double v, double w)
{
    int n, t, i, j, k, m;
    const double *x, *y, *z;
    double a[3], b[3], c[3], d[3], e[3];
    double zm, eps;

    if (q->Update == 0)
      ERR(("surface_update was not called"));
    eps = 1e-10;
    x = q->x;
    y = q->y;
    z = q->z;
    //    vec_ini(u, v, max(zm, w) + eps, /**/ e);
    //n = he_nt(he);
    for (t = m = 0; t < n; t++) {
      //he_tri_ijk(he, t, &i, &j, &k);
      //vec_get(i, x, y, z, a);
      //vec_get(j, x, y, z, b);
      //vec_get(k, x, y, z, c);
      //m += predicate_ray(d, e, a, b, c);
    }
    return m % 2;
}
