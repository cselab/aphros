// Created by Sergey Litvinov on 24.12.2019
// Copyright 2019 ETH Zurich

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>
#include <aphros/march/march.h>

enum { X, Y, Z };
static double r = 1.0;
static double lo = -1.1, hi = 1.1;

static double f(double x, double y, double z) {
  return x * x + y * y + z * z - r * r;
}

static void write(double x, double y, double z, double d, int n, double* tri) {
  int i, j, u, v, w;
  double a, b, c;
  static int J = 1;

  if (n == 0) return;
  for (i = j = 0; i < 3 * n; i++) {
    a = d * tri[j++] + x;
    b = d * tri[j++] + y;
    c = d * tri[j++] + z;
    printf("v %g %g %g\n", a, b, c);
  }
  for (i = 0; i < n; i++) {
    u = J++;
    v = J++;
    w = J++;
    printf("f %d %d %d\n", u, v, w);
  }
}

int main(int argc, char** argv) {
  int n, m, i, j, k, l;
  double x, y, z, d, u[8], xx[8], yy[8], zz[8], tri[3 * 3 * MARCH_NTRI], *o;
  char* me;

  m = 0;
  me = *argv;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
      case 'h':
        fprintf(stderr, "%s -n int > OBJ\n", me);
        exit(2);
        break;
      case 'n':
        argv++;
        if (*argv == NULL) {
          fprintf(stderr, "%s: -n needs an argument\n", me);
          exit(2);
        }
        m = atoi(*argv);
        break;
      default:
        fprintf(stderr, "%s: unknow option\n", me);
        exit(2);
    }
  if (m == 0) {
    fprintf(stderr, "%s: -n is not set\n", me);
    exit(2);
  }

  printf("# File type: ASCII OBJ\n");
  d = (hi - lo) / (m - 1);
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++)
      for (k = 0; k < m; k++) {
        x = lo + d * i;
        y = lo + d * j;
        z = lo + d * k;
        for (l = 0; l < 8; l++) {
          o = MARCH_O[l];
          xx[l] = x + d * o[X];
          yy[l] = y + d * o[Y];
          zz[l] = z + d * o[Z];
          u[l] = f(xx[l], yy[l], zz[l]);
        }
        march_cube(u, &n, tri);
        write(x, y, z, d, n, tri);
      }
}
