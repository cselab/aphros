// Created by Sergey Litvinov on 31.01.2021
// Copyright 2021 ETH Zurich

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bbox.h"
#include "err.h"
#include "inside.h"
#include "memory.h"
#include "predicate.h"

#define SIZE (9999)

int off_read(FILE*, int* status, int*, int**, int*, double**);
int ply_read(FILE*, int* status, int*, int**, int*, double**);
int stl_read(FILE*, int* status, int*, int**, int*, double**);

static double max_edg(const double* u, const double* v, const double* w);
static double sq(double);
static double edg(const double*, const double*);
static double tri_point_distance2(
    const double[3], const double[3], const double[3], const double p[3]);
static double max_dbl(double, double);

static char me[] = "inside";
enum { X, Y, Z };

struct Inside {
  int nt;
  struct Bbox* bbox;
  const int* tri;
  const double* ver;

  struct {
    double lo[2];
    double size;
    int* cap;
    int* n;
    int** data;
    int nx;
    int ny;
  } list;
};

int inside_fin(struct Inside* q) {
  int i;
  bbox_fin(q->bbox);
  for (i = 0; i < q->list.nx * q->list.ny; i++)
    FREE(q->list.data[i]);
  FREE(q->list.cap);
  FREE(q->list.data);
  FREE(q->list.n);
  FREE(q);
  return 0;
}

int inside_ini(int nt, const int* tri, const double* ver, struct Inside** pq) {
  const double* a;
  const double* b;
  const double* c;
  const double* hi;
  const double* lo;
  const double* r;
  double diameter;
  double size;
  int* cap;
  int** data;
  int dx;
  int dy;
  int i;
  int ix;
  int iy;
  int j;
  int jx;
  int jy;
  int k;
  int m;
  int* n;
  int nx;
  int ny;
  int t;
  struct Bbox* bbox;
  struct Inside* q;

  MALLOC(1, &q);
  bbox_ini(&bbox);
  predicate_ini();
  m = 0;
  for (i = 0; i < 3 * nt; i++)
    if (tri[i] > m) m = tri[i];
  bbox_update(bbox, m, ver);
  size = 0;
  for (t = 0; t < nt; t++) {
    i = tri[3 * t];
    j = tri[3 * t + 1];
    k = tri[3 * t + 2];
    a = &ver[3 * i];
    b = &ver[3 * j];
    c = &ver[3 * k];
    diameter = max_edg(a, b, c);
    if (diameter > size) size = diameter;
  }
  bbox_lo(bbox, &lo);
  bbox_hi(bbox, &hi);
  nx = (hi[X] - lo[X] + size) / size;
  ny = (hi[Y] - lo[Y] + size) / size;
  MALLOC(nx * ny, &cap);
  MALLOC(nx * ny, &n);
  MALLOC(nx * ny, &data);
  for (i = 0; i < nx * ny; i++) {
    cap[i] = 1;
    n[i] = 0;
    MALLOC(1, &data[i]);
  }

  for (t = 0; t < nt; t++) {
    i = tri[3 * t];
    r = &ver[3 * i];
    ix = (r[X] - lo[X]) / size;
    iy = (r[Y] - lo[Y]) / size;
    for (dx = -1; dx < 2; dx++)
      for (dy = -1; dy < 2; dy++) {
        jx = ix + dx;
        jy = iy + dy;
        j = jx + jy * nx;
        if (j < 0) continue;
        if (j >= nx * ny) continue;
        if (n[j] >= cap[j]) {
          cap[j] *= 2;
          REALLOC(cap[j], &data[j]);
        }
        data[j][n[j]++] = t;
      }
  }
  q->nt = nt;
  q->tri = tri;
  q->ver = ver;
  q->bbox = bbox;
  q->list.lo[X] = lo[X];
  q->list.lo[Y] = lo[Y];
  q->list.size = size;
  q->list.n = n;
  q->list.cap = cap;
  q->list.data = data;
  q->list.nx = nx;
  q->list.ny = ny;
  *pq = q;
  return 0;
}

int inside_inside_naive(struct Inside* q, const double r[3]) {
  int nt;
  int t;
  int i;
  int j;
  int k;
  int m;
  const int* tri;
  const double* ver;
  const double* a;
  const double* b;
  const double* c;
  double e[3];
  double zm;
  double eps;
  struct Bbox* bbox;

  eps = 1e-10;
  ver = q->ver;
  nt = q->nt;
  tri = q->tri;
  bbox = q->bbox;
  if (!bbox_inside(bbox, r)) {
    return 0;
  }
  zm = bbox_zhi(bbox);
  e[X] = r[X];
  e[Y] = r[Y];
  e[Z] = max_dbl(zm, r[Z]) + eps;
  for (t = m = 0; t < nt; t++) {
    i = *tri++;
    j = *tri++;
    k = *tri++;
    a = &ver[3 * i];
    b = &ver[3 * j];
    c = &ver[3 * k];
    m += predicate_ray(r, e, a, b, c);
  }
  return m % 2;
}

int inside_inside(struct Inside* q, const double r[3]) {
  const double* a;
  const double* b;
  const double* c;
  const double* lo;
  const double* ver;
  const int* tri;
  double e[3];
  double eps;
  double size;
  double zm;
  int** data;
  int i;
  int idx;
  int intersect;
  int item;
  int ix;
  int iy;
  int j;
  int k;
  int* n;
  int nx;
  int ny;
  int t;
  struct Bbox* bbox;

  eps = 1e-10;
  ver = q->ver;
  tri = q->tri;
  bbox = q->bbox;
  nx = q->list.nx;
  ny = q->list.ny;
  data = q->list.data;
  size = q->list.size;
  n = q->list.n;
  lo = q->list.lo;
  if (!bbox_inside(bbox, r)) {
    return 0;
  }
  zm = bbox_zhi(bbox);
  e[X] = r[X];
  e[Y] = r[Y];
  e[Z] = max_dbl(zm, r[Z]) + eps;
  ix = (r[X] - lo[X]) / size;
  iy = (r[Y] - lo[Y]) / size;
  idx = ix + iy * nx;
  if (idx < 0) idx = 0;
  if (idx >= nx * ny) idx = nx * ny - 1;
  intersect = 0;
  for (item = 0; item < n[idx]; item++) {
    t = data[idx][item];
    i = tri[3 * t];
    j = tri[3 * t + 1];
    k = tri[3 * t + 2];
    a = &ver[3 * i];
    b = &ver[3 * j];
    c = &ver[3 * k];
    intersect += predicate_ray(r, e, a, b, c);
  }
  return intersect % 2;
}

int inside_fwrite(struct Inside* q, FILE* file) {
  int nx;
  int ny;
  int ix;
  int iy;
  int idx;
  const int* n;
  nx = q->list.nx;
  ny = q->list.ny;
  n = q->list.n;
  for (ix = 0; ix < nx; ix++)
    for (iy = 0; iy < ny; iy++) {
      idx = ix + iy * nx;
      if (idx < 0) idx = 0;
      if (idx >= nx * ny) idx = nx * ny - 1;
      if (fprintf(file, "%d %d %d\n", ix, iy, n[idx]) < 0) return 1;
    }
  return 0;
}

int inside_info(struct Inside* q, struct InsideInfo* info) {
  int nx;
  int ny;
  int ix;
  int iy;
  int idx;
  int tri;
  int max_tri;
  int min_tri;
  const int* n;
  nx = info->nx = q->list.nx;
  ny = info->ny = q->list.ny;
  info->size = q->list.size;
  n = q->list.n;
  max_tri = 0;
  min_tri = INT_MAX;
  for (ix = 0; ix < nx; ix++)
    for (iy = 0; iy < ny; iy++) {
      idx = ix + iy * nx;
      if (idx < 0) idx = 0;
      if (idx >= nx * ny) idx = nx * ny - 1;
      tri = n[idx];
      if (tri > max_tri) max_tri = tri;
      if (tri < min_tri) min_tri = tri;
    }

  info->min_tri = min_tri;
  info->max_tri = max_tri;

  return 0;
}

double inside_distance(struct Inside* q, const double r[3]) {
  const double* a;
  const double* b;
  const double* c;
  const double* lo;
  const double* ver;
  const int* tri;
  double d;
  double mi;
  double size;
  int** data;
  int i;
  int idx;
  int item;
  int ix;
  int iy;
  int j;
  int k;
  int* n;
  int nx;
  int ny;
  int t;

  ver = q->ver;
  tri = q->tri;
  nx = q->list.nx;
  ny = q->list.ny;
  data = q->list.data;
  size = q->list.size;
  n = q->list.n;
  lo = q->list.lo;
  ix = (r[X] - lo[X]) / size;
  iy = (r[Y] - lo[Y]) / size;
  idx = ix + iy * nx;
  if (idx < 0) idx = 0;
  if (idx >= nx * ny) idx = nx * ny - 1;

  mi = DBL_MAX;
  for (item = 0; item < n[idx]; item++) {
    t = data[idx][item];
    i = tri[3 * t];
    j = tri[3 * t + 1];
    k = tri[3 * t + 2];
    a = &ver[3 * i];
    b = &ver[3 * j];
    c = &ver[3 * k];
    d = tri_point_distance2(a, b, c, r);
    if (d < mi) mi = d;
  }
  d = sqrt(mi);
  if (inside_inside(q, r))
    return -d;
  else
    return d;
}

double inside_distance_naive(struct Inside* q, const double r[3]) {
  const double* a;
  const double* b;
  const double* c;
  const double* ver;
  const int* tri;
  double d;
  double mi;
  int i;
  int j;
  int k;
  int t;
  int nt;

  ver = q->ver;
  tri = q->tri;
  nt = q->nt;
  mi = DBL_MAX;
  for (t = 0; t < nt; t++) {
    i = tri[3 * t];
    j = tri[3 * t + 1];
    k = tri[3 * t + 2];
    a = &ver[3 * i];
    b = &ver[3 * j];
    c = &ver[3 * k];
    d = tri_point_distance2(a, b, c, r);
    if (d < mi) mi = d;
  }
  d = sqrt(mi);
  if (inside_inside(q, r))
    return -d;
  else
    return d;
}

typedef int (*const ReadType)(
    FILE*, int* status, int* nt, int** tri, int* nv, double** ver);
static ReadType Read[] = {off_read, ply_read, stl_read};
int inside_mesh_read(
    const char* path, int* nt, int** tri, int* nv, double** ver) {
  int status;
  int err;
  FILE* file;
  long unsigned int i;

  for (i = 0; i < sizeof(Read) / sizeof(Read[0]); i++) {
    if ((file = fopen(path, "r")) == NULL) goto err;
    err = Read[i](file, &status, nt, tri, nv, ver);
    if (err != 0) goto err;
    fclose(file);
    if (status == 0) goto ok;
  }
err:
  return 1;
ok:
  return 0;
}

int inside_mesh_fin(int* tri, double* ver) {
  free(tri);
  free(ver);
  return 0;
}

int inside_box(struct Inside* q, double lo[3], double hi[3]) {
  const double* d;
  bbox_lo(q->bbox, &d);
  lo[X] = d[X];
  lo[Y] = d[Y];
  lo[Z] = d[Z];

  bbox_hi(q->bbox, &d);
  hi[X] = d[X];
  hi[Y] = d[Y];
  hi[Z] = d[Z];
  return 0;
}

static double sq(double x) {
  return x * x;
}
static double edg(const double* a, const double* b) {
  enum { X, Y, Z };
  return sqrt(sq(a[X] - b[X]) + sq(a[Y] - b[Y]) + sq(a[Z] - b[Z]));
}
double max_edg(const double* u, const double* v, const double* w) {
  double a;
  double b;
  double c;
  a = edg(v, u);
  b = edg(w, u);
  c = edg(v, w);
  return max_dbl(c, max_dbl(a, b));
}

static int vec_minus(const double a[3], const double b[3], /**/ double c[3]) {
  c[X] = a[X] - b[X];
  c[Y] = a[Y] - b[Y];
  c[Z] = a[Z] - b[Z];
  return 0;
}

static double vec_dot(const double a[3], const double b[3]) {
  return a[X] * b[X] + a[Y] * b[Y] + a[Z] * b[Z];
}

static double edg_sq(const double a[3], const double b[3]) {
  double u[3];

  vec_minus(b, a, u);
  return vec_dot(u, u);
}

static double edg_point_distance2(
    const double a[3], const double b[3], const double p[3]) {
  double t, s, x, y, z;

  s = edg_sq(a, b);
  if (s == 0) return edg_sq(p, a);
  t = ((b[X] - a[X]) * (p[X] - a[X]) + (b[Y] - a[Y]) * (p[Y] - a[Y]) +
       (b[Z] - a[Z]) * (p[Z] - a[Z])) /
      s;
  if (t > 1.0) return edg_sq(p, b);
  if (t < 0.0) return edg_sq(p, a);
  x = (1 - t) * a[X] + t * b[X] - p[X];
  y = (1 - t) * a[Y] + t * b[Y] - p[Y];
  z = (1 - t) * a[Z] + t * b[Z] - p[Z];
  return x * x + y * y + z * z;
}

static double tri_point_distance2(
    const double a[3], const double b[3], const double c[3],
    const double p[3]) {
  enum { X, Y, Z };
  double u[3], v[3], q[3];
  double A, B, C, D, E, det;
  double t1, t2;
  double x, y, z;
  double d1, d2;

  vec_minus(b, a, u);
  vec_minus(c, a, v);
  B = vec_dot(v, u);
  E = vec_dot(u, u);
  C = vec_dot(v, v);
  det = B * B - E * C;
  if (det == 0) {
    d1 = edg_point_distance2(a, b, p);
    d2 = edg_point_distance2(b, c, p);
    if (d1 < d2) return d1;
    return d2;
  }
  vec_minus(a, p, q);
  A = vec_dot(v, q);
  D = vec_dot(u, q);
  t1 = (D * C - A * B) / det;
  t2 = (A * E - D * B) / det;
  if (t1 < 0) return edg_point_distance2(a, c, p);
  if (t2 < 0) return edg_point_distance2(a, b, p);
  if (t1 + t2 > 1) return edg_point_distance2(b, c, p);
  x = q[X] + t1 * u[X] + t2 * v[X];
  y = q[Y] + t1 * u[Y] + t2 * v[Y];
  z = q[Z] + t1 * u[Z] + t2 * v[Z];
  return x * x + y * y + z * z;
}

static double max_dbl(double a, double b) {
  return a > b ? a : b;
}
