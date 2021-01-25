#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <csv.h>
#include <table.h>

enum { N = 999 };
static const char me[] = "csv2sph";

#include "ico.inc"
#include "util.h"

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
  char output[N];
  char* Prefix;
  double* r;
  double* x;
  double* y;
  double* z;
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
    for (i = 0; i < nr; i++)
      for (j = 0; j < ico_nv; j++)
        fprintf(
            file, "%.16g %.16g %.16g\n", x[i] + r[i] * ico_x[j],
            y[i] + r[i] * ico_y[j], z[i] + r[i] * ico_z[j]);
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
