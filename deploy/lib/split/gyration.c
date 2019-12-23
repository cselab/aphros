#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <chmath.h>
#include <csv.h>
#include <table.h>

enum { N = 999 };
static const char me[] = "gyration";

#include "util.h"

#define GET(f, r)                                                   \
  if ((*r = csv_field(csv, f)) == NULL) {                           \
    fprintf(stderr, "%s: no field '%s' in '%s'\n", me, (f), *argv); \
    exit(2);                                                        \
  }

#define ADD(f, r)                                                        \
  do {                                                                   \
    if (csv_add(csv, f) != 0) {                                          \
      fprintf(stderr, "%s: fail to add '%s' to '%s'\n", me, (f), *argv); \
      exit(2);                                                           \
    }                                                                    \
    GET(f, r);                                                           \
  } while (0)

#define USED(x) \
  if (x)        \
    ;           \
  else {        \
  }
#define MALLOC(n, p)                                      \
  do {                                                    \
    *(p) = malloc((n) * sizeof(**(p)));                   \
    if (*(p) == NULL) {                                   \
      fprintf(stderr, "%s: alloc failed, n = %d", me, n); \
      exit(2);                                            \
    }                                                     \
  } while (0)
#define REALLOC(n, p)                                       \
  do {                                                      \
    *(p) = realloc(*(p), (n) * sizeof(**(p)));              \
    if (*(p) == NULL) {                                     \
      fprintf(stderr, "%s: realloc failed, n = %d", me, n); \
      exit(2);                                              \
    }                                                       \
  } while (0)

static void usg() {
  fprintf(stderr, "%s -p prefix [csv..]\n", me);
  exit(1);
}

static int gyration(
    double x, double, double, double xx, double, double, double, double, double,
    double*, double*, double*);

int main(int argc, char** argv) {
  enum { X, Y, Z };
  char output[N];
  char* Prefix;
  double* asphericity;
  double* acylindricity;
  double* lx;
  double* ly;
  double* lz;
  double* rg;
  double* x;
  double* xx;
  double* xy;
  double* xz;
  double* y;
  double* yy;
  double* yz;
  double* z;
  double* zz;
  FILE* file;
  int i;
  int nr;
  struct CSV* csv;

  USED(argc);
  Prefix = NULL;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
      case 'h':
        usg();
        break;
      case 'p':
        argv++;
        if ((Prefix = *argv) == NULL) {
          fprintf(stderr, "%s: -p needs an argument\n", me);
          exit(2);
        }
        break;
      default:
        fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
        exit(1);
    }
  if (Prefix == NULL) {
    fprintf(stderr, "%s: prefix (-p) is not given\n", me);
    exit(1);
  }
  if (*argv == NULL) {
    fprintf(stderr, "%s: csv file is not given\n", me);
    exit(1);
  }

  for (; *argv != NULL; argv++) {
    if ((file = fopen(*argv, "r")) == NULL) {
      fprintf(stderr, "%s: fail to open '%s'\n", me, *argv);
      exit(1);
    }
    if ((csv = csv_read(file)) == NULL) {
      fprintf(stderr, "%s: not a cvs file '%s'\n", me, *argv);
      exit(1);
    }
    fclose(file);

    GET("x", &x);
    GET("y", &y);
    GET("z", &z);
    GET("xx", &xx);
    GET("xy", &xy);
    GET("xz", &xz);
    GET("yy", &yy);
    GET("yz", &yz);
    GET("zz", &zz);
    ADD("lx", &lx);
    ADD("ly", &ly);
    ADD("lz", &lz);
    ADD("asphericity", &asphericity);
    ADD("rg", &rg);
    ADD("acylindricity", &acylindricity);
    nr = csv_nr(csv);
    double a, b, c;

    for (i = 0; i < nr; i++) {
      gyration(
          x[i], y[i], z[i], xx[i], xy[i], xz[i], yy[i], yz[i], zz[i], &a, &b,
          &c);
      lx[i] = a;
      ly[i] = b;
      lz[i] = c;
      rg[i] = sqrt(a + b + c);
      asphericity[i] = c - (a + b) / 2;
      acylindricity[i] = b - a;
    }
    if (util_name(Prefix, *argv, output) != 0) {
      fprintf(stderr, "%s: util_name failed\n", me);
      exit(2);
    }
    if ((file = fopen(output, "w")) == NULL) {
      fprintf(stderr, "%s: fail to write to '%s'\n", me, output);
      exit(2);
    }
    if (csv_write(csv, file) != 0) {
      fprintf(stderr, "%s: fail to wrote to '%s'\n", me, *argv);
      exit(2);
    }
    fclose(file);
  }
}

static int gyration(
    double x, double y, double z, double xx, double xy, double xz, double yy,
    double yz, double zz, double* a, double* b, double* c) {
  enum { X, Y, Z };
  enum { XX, XY, XZ, YY, YZ, ZZ };
  double i[6];
  double val[3];

  i[XX] = xx - x * x;
  i[XY] = xy - x * y;
  i[XZ] = xz - x * z;
  i[YY] = yy - y * y;
  i[YZ] = yz - y * z;
  i[ZZ] = zz - z * z;
  if (chmath_eig_values(i, val) != 0) {
    fprintf(stderr, "%s: math_eig_values failed\n", me);
    exit(2);
  }
  *a = val[X];
  *b = val[Y];
  *c = val[Z];
  return 0;
}
