#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inside.h>

static const char* me = "plane";
static double eps = 1e-10;
#define USED(x) \
  if (x)        \
    ;           \
  else {        \
  }

static void usg(void) {
  fprintf(stderr, "%s mesh > off\n", me);
  exit(1);
}

int main(int argc, const char** argv) {
  USED(argc);
  int nt;
  int nv;
  int* tri;
  double* ver;

  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
      case 'h':
        usg();
        break;
      default:
        fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
        exit(2);
    }
  if (argv[0] == NULL) {
    fprintf(stderr, "%s: mesh file is missing\n", me);
    exit(2);
  }
  if (inside_mesh_read(argv[0], &nt, &tri, &nv, &ver) != 0) {
    fprintf(stderr, "%s: fail to read mesh '%s'\n", me, argv[0]);
    exit(2);
  }
  inside_mesh_fin(tri, ver);
}

static double plane_point(const double n[3], double alpha, const double p[3])
{
  enum {X, Y, Z};
  return n[X]*p[X] + n[Y]*p[Y] + n[Z]*p[Z] + alpha;
}

static int plain_edg(const double n[3], double alpha, const double a[3], const double b[3], /**/ double p[3])
{
  enum {X, Y, Z};  
  double x;
  double y;
  double t;
  x = plane_point(n, alpha, a);
  if (fabs(x) < eps) {
    p[X] = a[X];
    p[Y] = a[Y];
    p[Z] = a[Z];
    return 1;
  }
  y = plane_point(n, alpha, b);
  if (fabs(y) < eps) {
    p[X] = b[X];
    p[Y] = b[Y];
    p[Z] = b[Z];
    return 1;
  }  
  if (x * y > 0)
    return 0;
  t = x / (x - y);
  p[X] = a[X] + t * (b[X] - a[X]);
  p[Y] = a[Y] + t * (b[Y] - a[Y]);
  p[Y] = a[Z] + t * (b[Z] - a[Z]);
  return 1;
}

static int plain_tri(const double n[3], double alpha, const double a[3], const double b[3], const double c[3], double *ver) {
  int cnt;
  cnt = 0;
  if (plain_edg(n, alpha, a, b, ver)) {
    cnt ++;
    ver += 3;
  }
  if (plain_edg(n, alpha, b, c, ver)) {
    cnt ++;
    ver += 3;
  }
  if (plain_edg(n, alpha, c, a, ver)) {
    cnt ++;
    ver += 3;
  }  
  return cnt;
}
