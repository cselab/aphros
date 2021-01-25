#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include <march.h>

#define SIZE(x) (sizeof(x) / sizeof(*(x)))

enum { X, Y, Z };
static double r = 0.25;
static double lo = -0.5, hi = 0.5;
static int m = 10;
static double av(double a, double b, double o) {
  return a + (b - a) * o;
}

static double sq(double x) {
  return x * x;
}

static double f(double x, double y, double z) {
  double ans;

  ans = sq(x) + sq(2 * y) + sq(3 * z) - sq(r);
  return ans;
}

static double O[][3] = {
    {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
    {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1},
};

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
  double tri[3 * 3 * MARCH_NTRI];
  double offset[3 * MARCH_NTRI];
  int p[3 * MARCH_NTRI], q[3 * MARCH_NTRI];
  int n, i, j, k, l, u;
  int a, b, w;
  double x, y, z, d, pos;
  double cube[8], xx[8], yy[8], zz[8];
  double* o;
  int stat[MARCH_NTRI] = {0};

  while (argv[1] != NULL && argv++[1][0] == '-')
    switch (argv[0][1]) {
      case 'n':
        if (argv[1] == NULL) {
          fprintf(stderr, "%s: -n needs an argument\n", argv[0]);
          exit(2);
        }
        m = atoi(argv++[1]);
        break;
      default:
        fprintf(stderr, "%s: unknow option\n", argv[0]);
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
          o = O[l];
          xx[l] = x + d * o[X];
          yy[l] = y + d * o[Y];
          zz[l] = z + d * o[Z];
          cube[l] = f(xx[l], yy[l], zz[l]);
        }
        /* march_cube(cube, &n, tri);
           march_cube_location(p, q, offset); */
        march_cube_location(cube, &n, tri, p, q, offset);
        for (u = 0; u < n; u++) {
          pos = av(O[p[u]][X], O[q[u]][X], offset[u]);
          assert(fabs(pos - tri[3 * u + X] < 1e-12));

          pos = av(O[p[u]][Y], O[q[u]][Y], offset[u]);
          assert(fabs(pos - tri[3 * u + Y] < 1e-12));

          pos = av(O[p[u]][Z], O[q[u]][Z], offset[u]);
          assert(fabs(pos - tri[3 * u + Z] < 1e-12));
        }
        stat[n]++;
        write(x, y, z, d, n, tri);
      }
  for (i = 0; i < SIZE(stat); i++)
    if (stat[i] > 0) fprintf(stderr, "%2d %7d\n", i, stat[i]);
}
