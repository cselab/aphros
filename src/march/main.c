// Created by Petr Karnakov on 31.01.2021
// Copyright 2021 ETH Zurich

#include <stdio.h>

#include "march.h"
#include "table.inc"

enum { X, Y, Z };
struct Ver {
  double a;
  int x, y;
};
struct March {
  struct Ver cube_ver[3 * MARCH_NTRI];
  int cube_n;
};

static double offset(double, double);
static int map(int);
const int MARCH_O[][3] = {
    {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
    {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1},
};
static int cube(struct March* q, double cube[8], int* pn, double* tri) {
  double a;
  int c, i, j, idx, flag, x, y;
  double *v, *o, *dir;
  int n, k, m;
  struct {
    int x, y;
    double a;
    double v[3];
  } ve[12];

  idx = 0;
  for (i = 0; i < 8; i++) {
    if (cube[i] <= 0) idx |= 1 << i;
  }
  if ((flag = CubeEdgeFlags[idx]) == 0) {
    *pn = q->cube_n = 0;
    return 0;
  }
  for (i = 0; i < 12; i++) {
    if (flag & (1 << i)) {
      o = Offset[Connection[i][0]];
      dir = Direction[i];
      ve[i].x = x = Connection[i][0];
      ve[i].y = y = Connection[i][1];
      ve[i].a = a = offset(cube[x], cube[y]);
      ve[i].v[X] = o[X] + a * dir[X];
      ve[i].v[Y] = o[Y] + a * dir[Y];
      ve[i].v[Z] = o[Z] + a * dir[Z];
    }
  }

  n = k = m = 0;
  for (i = 0; i < 5; i++) {
    if (TriangleConnectionTable[idx][3 * i] < 0) break;
    for (c = 0; c < 3; c++) {
      j = TriangleConnectionTable[idx][3 * i + c];
      v = ve[j].v;
      tri[k++] = v[X];
      tri[k++] = v[Y];
      tri[k++] = v[Z];
      q->cube_ver[m].x = ve[j].x;
      q->cube_ver[m].y = ve[j].y;
      q->cube_ver[m].a = ve[j].a;
      m++;
    }
    n++;
  }
  *pn = q->cube_n = n;
  return 0;
}

static double offset(double a, double b) {
  double d;

  d = a - b;
  if (d == 0.0) return 0.5;
  return a / d;
}

static void swap(double* u, int i, int j) {
  double t;

  t = u[i];
  u[i] = u[j];
  u[j] = t;
}

static int march(struct March* q, double u[8], int* pn, double* tri) {
  int s;

  swap(u, 2, 3);
  swap(u, 6, 7);
  s = cube(q, u, pn, tri);
  swap(u, 2, 3);
  swap(u, 6, 7);
  return s;
}

int march_cube(double u[8], int* n, double* tri) {
  struct March q;
  return march(&q, u, n, tri);
}

static int cube_location(struct March* q, int* x, int* y, double* a) {
  int i;

  for (i = 0; i < 3 * q->cube_n; i++) {
    x[i] = map(q->cube_ver[i].x);
    y[i] = map(q->cube_ver[i].y);
    a[i] = q->cube_ver[i].a;
  }
  return 0;
}

int march_cube_location(
    double u[8], /**/ int* ntri, double* tri, int* x, int* y, double* a) {
  int status;
  struct March q;

  status = march(&q, u, ntri, tri);
  if (status != 0) return status;
  status = cube_location(&q, x, y, a);
  return status;
}

static int map(int i) {
  int m[8] = {0, 1, 3, 2, 4, 5, 7, 6};
  return m[i];
}
