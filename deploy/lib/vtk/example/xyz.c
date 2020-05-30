#include <stdio.h>
#include <stdlib.h>

#include <vtk.h>

#define USED(x) \
  if (x)        \
    ;           \
  else {        \
  }
static char me[] = "vtk/xyz";
static void usg() {
  fprintf(stderr, "%s [-b] < vtk\n", me);
  exit(1);
}

static double array_min(int, double[]);
static double array_max(int, double[]);

int main(int argc, char** argv) {
  USED(argc);
  struct VTK* vtk;
  int nv, i, Box;
  double xl, xh, yl, yh, zl, zh;

  Box = 0;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
      case 'h':
        usg();
        break;
      case 'b':
        argv++;
        Box = 1;
        break;
      default:
        fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
        exit(1);
    }
  vtk = vtk_read(stdin);
  if (vtk == NULL) {
    fprintf(stderr, "%s: fail to read\n", me);
    exit(2);
  }
  nv = vtk_nv(vtk);

  if (Box) {
    xh = array_max(nv, vtk->x);
    xl = array_min(nv, vtk->x);
    yh = array_max(nv, vtk->y);
    yl = array_min(nv, vtk->y);
    zh = array_max(nv, vtk->z);
    zl = array_min(nv, vtk->z);
    printf("%g %g %g %g %g %g\n", xl, xh, yl, yh, zl, zh);
  } else {
    for (i = 0; i < nv; i++)
      printf("%.16g %.16g %.16g\n", vtk->x[i], vtk->y[i], vtk->z[i]);
  }

  vtk_fin(vtk);
}

static double array_max(int n, double a[]) {
  int i;
  double m;

  m = a[0];
  for (i = 0; i < n; i++) {
    if (a[i] > m) m = a[i];
  }
  return m;
}

static double array_min(int n, double a[]) {
  int i;
  double m;

  m = a[0];
  for (i = 0; i < n; i++) {
    if (a[i] < m) m = a[i];
  }
  return m;
}
