#include <math.h>
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

int off_read(FILE *, int * status, int *, int **, int *, double **);
int ply_read(FILE *, int * status, int *, int **, int *, double **);

static int circumradius(const double * u, const double * v, const double * w, double *r);
static double sq(double);
static double edg(const double *, const double *);

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
  int t;
  int i;
  int j;
  int k;
  int m;
  const double *a;
  const double *b;
  const double *c;
  struct Inside *q;
  double radii;
  double diameter;
  double max_diameter;

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

  max_diameter = 0;
  for (t = 0; t < n; t++) {
    i = tri[3*t];
    j = tri[3*t+1];
    k = tri[3*t+2];
    a = &ver[3*i];
    b = &ver[3*j];
    c = &ver[3*k];
    if (circumradius(a, b, c, &radii) != 0)
      return 1;
    diameter = 2 * radii;
    if (diameter > max_diameter)
      max_diameter = diameter;
  }
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
  double e[3];
  double zm;
  double eps;
  struct Bbox * bbox;

  eps = 1e-10;
  ver = q->ver;
  n = q->n;
  tri = q->tri;
  bbox = q->bbox;
  if (!bbox_inside(bbox, r)) {
    return 0;
  }
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
    m += predicate_ray(r, e, a, b, c);
  }
  return m % 2;
}

typedef int (*const ReadType)(FILE *, int * status, int * nt, int ** tri, int * nv, double ** ver);
static const ReadType Read[] = {off_read, ply_read};
 int
inside_mesh_read(const char *path, int * nt, int ** tri, int * nv, double ** ver)
{
  int status;
  int err;
  int state;
  ReadType read;
  FILE *file;
  long unsigned int i;

  for (i = 0; i < sizeof(Read)/sizeof(Read[0]); i++) {
    if ((file = fopen(path, "r")) == NULL)
      goto err;
    err = Read[i](file, &status, nt, tri, nv, ver);
    if (err != 0)
      goto err;
    fclose(file);
    if (status == 0)
      goto ok;
  }
 err:
  return 1;
 ok:
  return 0;
}

int
inside_mesh_fin(int * tri, double * ver)
{
  free(tri);
  free(ver);
  return 0;
}

static double
sq(double x)
{
  return x*x;
}
static double
edg(const double *a, const double *b)
{
  enum {X, Y, Z};
  return sqrt(sq(a[X] - b[X]) + sq(a[Y] - b[Y]) + sq(a[Z] - b[Z]));
}
static int
circumradius(const double * u, const double * v, const double * w, double *r)
{
  double a;
  double b;
  double c;
  double s;
  double num;
  double den;
  a = edg(v, u);
  b = edg(w, u);
  c = edg(v, w);
  s = (a + b + c)/2;
  num = a*b*c;
  den = s*(s - a)*(s - b)*(s - c);
  if (den <= 0)
    return 1;
  *r = num/(4*sqrt(den));
  return 0;
}
