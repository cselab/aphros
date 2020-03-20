#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "bbox.h"
#include "err.h"
#include "memory.h"
#include "predicate.h"
#include "inside.h"

static char me[] = "inside";
static int get(int, const double[], const double[], const double[], /**/ double a[3]);
enum {
      X, Y, Z
};

struct Inside {
  int n;
  struct Bbox *bbox;
  const int *tri;
  const double *x;
  const double *y;
  const double *z;
  int Update;
};

int
inside_ini(double lo[2], double hi[2], double size, struct Inside ** pq)
{
  struct Inside *q;

  MALLOC(1, &q);
  bbox_ini(&q->bbox);
  predicate_ini();
  q->Update = 0;
  *pq = q;
  return 0;
}

int
inside_fin(struct Inside * q)
{
  bbox_fin(q->bbox);  
  FREE(q);
  return 0;
}

int
inside_update(struct Inside * q, int n, const int * tri, const double * x, const double * y,
	      const double * z)
{
  q->n = n;
  q->tri = tri;
  q->x = x;
  q->y = y;
  q->z = z;
  bbox_update(q->bbox, n, x, y, z);  
  q->Update = 1;
  return 0;
}

#define max(a, b) ( (a) > (b) ? (a) : (b) )
int
inside_inside(struct Inside * q, double u, double v, double w)
{
  int n;
  int t;
  int i;
  int j;
  int k;
  int m;
  const int *tri;
  const double *x;
  const double *y;
  const double *z;
  double a[3];
  double b[3];
  double c[3];
  double d[3];
  double e[3];
  double zm;
  double eps;
  struct Bbox * bbox;

  if (q->Update == 0)
    ERR(("inside_update was not called"));
  eps = 1e-10;
  x = q->x;
  y = q->y;
  z = q->z;
  n = q->n;
  tri = q->tri;
  bbox = q->bbox;
  if (!bbox_inside(bbox, u, v, w))
    return 1;
  zm = bbox_zhi(bbox);
  e[X] = u;
  e[Y] = v;
  e[Z] = max(zm, w) + eps;
  for (t = m = 0; t < n; t++) {
    i = *tri++;
    j = *tri++;
    k = *tri++;
    get(i, x, y, z, a);
    get(j, x, y, z, b);
    get(k, x, y, z, c);
    m += predicate_ray(d, e, a, b, c);
  }
  return m % 2;
}

static int
get(int i, const double x[], const double y[], const double z[], /**/ double a[3])
{
  a[X] = x[i];
  a[Y] = y[i];
  a[Z] = z[i];
  return 0;
}
