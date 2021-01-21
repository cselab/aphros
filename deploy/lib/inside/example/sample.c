#include <stdio.h>
#include <stdlib.h>
#include <inside.h>

static const char* me = "sample";
enum { seed = 12345 };
#define USED(x) \
  if (x)        \
    ;           \
  else {        \
  }
static void usg(void) {
  fprintf(stderr, "%s [-i] -p float float float < off\n", me);
  exit(1);
}
static double rnd(double, double);

int main(int argc, const char** argv) {
  USED(argc);
  enum { X, Y, Z };
  double a;
  double b;
  double hi[3];
  double lo[3];
  double r[3];
  double size;
  double* ver;
  int i;
  int Invert;
  int Nflag;
  int n;
  int nt;
  int nv;
  int tmp;
  int* tri;
  struct InsideInfo info;
  struct Inside* inside;

  Invert = 0;
  Nflag = 0;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
      case 'h':
        usg();
        break;
      case 'i':
        Invert = 1;
        break;
      case 'n':
        Nflag = 1;
        argv++;
        if (*argv == NULL) {
          fprintf(stderr, "%s: -n needs an argument\n", me);
          return 2;
        }
        n = atoi(*argv);
        break;
      default:
        fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
        return 2;
    }
  if (Nflag == 0) {
    fprintf(stderr, "%s: -n is not set\n", me);
    return 2;
  }
  if (argv[0] == NULL) {
    fprintf(stderr, "%s: mesh file is missing\n", me);
    return 2;
  }
  if (inside_mesh_read(argv[0], &nt, &tri, &nv, &ver) != 0) {
    fprintf(stderr, "%s: fail to read mesh '%s'\n", me, argv[0]);
    return 2;
  }
  if (Invert)
    for (i = 0; i < nt; i++) {
      tmp = tri[3 * i];
      tri[3 * i] = tri[3 * i + 1];
      tri[3 * i + 1] = tmp;
    }
  inside_ini(nt, tri, ver, &inside);
  inside_box(inside, lo, hi);
  inside_info(inside, &info);
  size = info.size;
  fprintf(stderr, "size: %g\n", size);
  srand(seed);
  for (i = 0; i < n; i++) {
    r[X] = rnd(lo[X], hi[X]);
    r[Y] = rnd(lo[Y], hi[Y]);
    r[Z] = rnd(lo[Z], hi[Z]);
    a = inside_distance(inside, r);
    b = inside_distance_naive(inside, r);
    if (-size < b && b < size && a != b)
      printf("%g %g %g %g %g\n", r[X], r[Y], r[Z], a, b);
  }

  inside_fin(inside);
  inside_mesh_fin(tri, ver);
}

static double rnd(double a, double b) {
  double s;
  s = (double)rand() / (RAND_MAX + 1.0);
  return s * (b - a) + a;
}
