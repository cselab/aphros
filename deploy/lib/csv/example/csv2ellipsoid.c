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
static const char me[] = "csv2ellipsoid";

#include "ico.inc"
#include "util.h"

struct Transform {
  double r[3];
  double a[3];
  double b[3];
  double c[3];
  double scale[3];
  double rg;
};

static int transform_ini(
    double x, double, double, double xx, double, double, double, double, double,
    struct Transform*);
static int transform_apply(
    struct Transform*, double x, double, double, double[3]);
static int transform_rg(struct Transform*, double*);

#define GET(f, r)                                                   \
  if ((*r = csv_field(csv, f)) == NULL) {                           \
    fprintf(stderr, "%s: no field '%s' in '%s'\n", me, (f), *argv); \
    exit(2);                                                        \
  }

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

int main(int argc, char** argv) {
  enum { X, Y, Z };
  char output[N];
  char* Prefix;
  double* r;
  double* x;
  double* y;
  double* z;
  double* xx;
  double* xy;
  double* xz;
  double* yy;
  double* yz;
  double* zz;
  double u[3];
  struct Transform transform;
  FILE* file;
  int i;
  int j;
  int k;
  int nf;
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
    GET("r", &r);
    nf = csv_nf(csv);
    nr = csv_nr(csv);

    if (util_name(Prefix, *argv, output) != 0) {
      fprintf(stderr, "%s: util_name failed\n", me);
      exit(2);
    }
    if ((file = fopen(output, "w")) == NULL) {
      fprintf(stderr, "%s: fail to write to '%s'\n", me, output);
      exit(2);
    }
    if (fputs("# vtk DataFile Version 2.0\n", file) == EOF) {
      fprintf(stderr, "%s: fail to write\n", me);
      return 1;
    }
    fprintf(file, "%s\n", me);
    fputs("ASCII\n", file);
    fputs("DATASET POLYDATA\n", file);
    fprintf(file, "POINTS %d double\n", nr * ico_nv);
    for (i = 0; i < nr; i++) {
      transform_ini(
          x[i], y[i], z[i], xx[i], xy[i], xz[i], yy[i], yz[i], zz[i],
          &transform);
      for (j = 0; j < ico_nv; j++) {
        transform_apply(&transform, ico_x[j], ico_y[j], ico_z[j], u);
        fprintf(file, "%.16g %.16g %.16g\n", u[X], u[Y], u[Z]);
      }
    }
    fprintf(file, "POLYGONS %d %d\n", nr * ico_nt, 4 * nr * ico_nt);
    for (i = 0; i < nr; i++)
      for (j = 0; j < ico_nt; j++)
        fprintf(
            file, "3 %d %d %d\n", ico_nv * i + ico_t0[j],
            ico_nv * i + ico_t1[j], ico_nv * i + ico_t2[j]);

    fprintf(file, "CELL_DATA %d\n", nr * ico_nt);
    for (i = 0; i < nf; i++) {
      fprintf(file, "SCALARS %s double 1\n", csv->name[i]);
      fprintf(file, "LOOKUP_TABLE default\n");
      for (j = 0; j < nr; j++)
        for (k = 0; k < ico_nt; k++)
          fprintf(file, "%.16g\n", csv->data[i][j]);
    }
    fclose(file);
    csv_fin(csv);
  }
}

int transform_ini(
    double x, double y, double z, double xx, double xy, double xz, double yy,
    double yz, double zz, struct Transform* t) {
  enum { X, Y, Z };
  enum { XX, XY, XZ, YY, YZ, ZZ };
  double i[6];
  double rg;
  double* scale;

  scale = t->scale;
  t->r[X] = x;
  t->r[Y] = y;
  t->r[Z] = z;
  i[XX] = xx - x * x;
  i[XY] = xy - x * y;
  i[XZ] = xz - x * z;
  i[YY] = yy - y * y;
  i[YZ] = yz - y * z;
  i[ZZ] = zz - z * z;
  rg = i[XX] + i[YY] + i[ZZ];
  rg = sqrt(5 * rg / 3);
  if (chmath_eig_vectors(i, t->a, t->b, t->c) != 0) {
    fprintf(stderr, "%s: math_eig_vectors failed\n", me);
    exit(2);
  }
  if (chmath_eig_values(i, t->scale) != 0) {
    fprintf(stderr, "%s: math_eig_values failed\n", me);
    exit(2);
  }
  t->rg = rg;
  scale[X] = sqrt(fabs(5 * scale[X]));
  scale[Y] = sqrt(fabs(5 * scale[Y]));
  scale[Z] = sqrt(fabs(5 * scale[Z]));
  return 0;
}

int transform_apply(
    struct Transform* t, double x, double y, double z, double v[3]) {
  enum { X, Y, Z };
  const double* a;
  const double* b;
  const double* c;
  const double* r;
  const double* scale;

  a = t->a;
  b = t->b;
  c = t->c;
  r = t->r;
  scale = t->scale;
  x *= scale[X];
  y *= scale[Y];
  z *= scale[Z];
  v[X] = x * a[X] + y * b[X] + z * c[X];
  v[Y] = x * a[Y] + y * b[Y] + z * c[Y];
  v[Z] = x * a[Z] + y * b[Z] + z * c[Z];

  v[X] += r[X];
  v[Y] += r[Y];
  v[Z] += r[Z];

  return 0;
}

static int transform_rg(struct Transform* t, double* p) {
  *p = t->rg;
  return 0;
}
