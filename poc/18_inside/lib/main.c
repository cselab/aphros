#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bbox.h"
#include "err.h"
#include "memory.h"
#include "predicate.h"
#include "inside.h"

#define SIZE (9999)

static char me[] = "inside";
enum {
      X, Y, Z
};

struct Inside {
  int n;
  struct Bbox *bbox;
  const int *tri;
  const double *ver;
};

int
inside_fin(struct Inside * q)
{
  bbox_fin(q->bbox);  
  FREE(q);
  return 0;
}

int
inside_ini(int n, const int * tri, const double * ver, struct Inside ** pq)
{
  int i;
  int m;
  struct Inside *q;
  
  MALLOC(1, &q);
  bbox_ini(&q->bbox);
  predicate_ini();
  q->n = n;
  q->tri = tri;
  q->ver = ver;
  m = 0;
  for (i = 0; i < 3*n; i++)
    if (tri[i] > m)
      m = tri[i];
  bbox_update(q->bbox, m, ver);
  *pq = q;
  return 0;
}

#define max(a, b) ( (a) > (b) ? (a) : (b) )
int
inside_inside(struct Inside * q, const double r[3])
{
  int n;
  int t;
  int i;
  int j;
  int k;
  int m;
  const int *tri;
  const double *ver;
  const double *a;
  const double *b;
  const double *c;  
  double d[3];
  double e[3];
  double zm;
  double eps;
  struct Bbox * bbox;

  eps = 1e-10;
  ver = q->ver;
  n = q->n;
  tri = q->tri;
  bbox = q->bbox;
  if (!bbox_inside(bbox, r))
    return 1;
  zm = bbox_zhi(bbox);
  e[X] = r[X];
  e[Y] = r[Y];
  e[Z] = max(zm, r[Z]) + eps;
  for (t = m = 0; t < n; t++) {
    i = *tri++;
    j = *tri++;
    k = *tri++;
    a = &ver[3*i];
    b = &ver[3*j];
    c = &ver[3*k];
    m += predicate_ray(d, e, a, b, c);
  }
  return m % 2;
}

