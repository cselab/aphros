#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <csv.h>
#include <table.h>
#include <vtk.h>

enum { N = 999 };
static char me[] = "vtk/color";

#include "util.h"

#define USED(x) \
  if (x)        \
    ;           \
  else {        \
  }
static double pi = 3.141592653589793;

static void usg() {
  fprintf(stderr, "%s -p prefix -k key -f field [csv ..] -- [vtk ..]\n", me);
  exit(1);
}

int main(int argc, char** argv) {
  int status, nt, nr, i, j, key;
  float* cl;
  double *cl_csv, *vf, *rad;
  struct VTK* vtk;
  struct CSV* csv;
  struct Table* table;
  char *path, *Prefix, *Key, *Volume, **Vtk, **Csv;
  char output[N];
  FILE* f;

  USED(argc);
  Prefix = Key = Volume = NULL;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
      case 'h':
        usg();
        break;
      case 'k':
        argv++;
        if ((Key = *argv) == NULL) {
          fprintf(stderr, "%s: -k needs an argument\n", me);
          exit(2);
        }
        break;
      case 'f':
        argv++;
        if ((Volume = *argv) == NULL) {
          fprintf(stderr, "%s: -f needs an argument\n", me);
          exit(2);
        }
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
  if (Key == NULL) {
    fprintf(stderr, "%s: key (-k) is not given\n", me);
    exit(1);
  }
  if (Volume == NULL) {
    fprintf(stderr, "%s: field (-f) is not given\n", me);
    exit(1);
  }
  if (Prefix == NULL) {
    fprintf(stderr, "%s: prefix (-p) is not given\n", me);
    exit(1);
  }
  Csv = argv;
  for (;;) {
    if (*argv == NULL) {
      fprintf(stderr, "%s: missing '--' in arguments\n", me);
      exit(2);
    }
    if (util_eq(*argv++, "--")) break;
  }
  Vtk = argv;
  for (;;) {
    if (*Vtk == NULL || *Csv == NULL) break;
    path = *Csv++;
    util_name(Prefix, path, output);
    f = fopen(path, "r");
    if (f == NULL) {
      fprintf(stderr, "%s: cannot read '%s'\n", me, path);
      exit(2);
    }
    csv = csv_read(f);
    if (csv == NULL) {
      fprintf(stderr, "%s: fail to parse csv file '%s'\n", me, path);
      exit(2);
    }
    fclose(f);
    path = *Vtk++;
    f = fopen(path, "r");
    if (f == NULL) {
      fprintf(stderr, "%s: cannot read '%s'\n", me, path);
      exit(2);
    }
    vtk = vtk_read(f);
    if (vtk == NULL) {
      fprintf(stderr, "%s: fail to parse csv file '%s'\n", me, path);
      exit(2);
    }
    fclose(f);

    nt = vtk_nt(vtk);
    cl = vtk_data(vtk, Key);
    if (cl == NULL) {
      fprintf(stderr, "%s: no field '%s' in vtk file\n", me, Key);
      exit(2);
    }
    cl_csv = csv_field(csv, Key);
    if (cl_csv == NULL) {
      fprintf(stderr, "%s: no field '%s' in csv file\n", me, Key);
      exit(2);
    }
    vf = csv_field(csv, Volume);
    if (vf == NULL) {
      fprintf(stderr, "%s: no field '%s' in csv file\n", me, Volume);
      exit(2);
    }
    nr = csv_nr(csv);
    table = table_ini(100);
    for (i = 0; i < nr; i++) {
      key = (int)cl_csv[i];
      table_put(table, key, i);
    }
    vtk_add(vtk, Volume, VTK_CELL, VTK_DOUBLE);
    rad = vtk_data(vtk, Volume);
    nt = vtk_nt(vtk);
    for (i = 0; i < nt; i++) {
      key = (int)cl[i];
      status = table_get(table, key, &j);
      if (status != TABLE_EMPY)
        rad[i] = vf[j];
      else
        rad[i] = 0;
    }
    if ((f = fopen(output, "w")) == NULL) {
      fprintf(stderr, "%s: fail to write to '%s'\n", me, output);
      exit(2);
    }
    vtk_write(vtk, f);
    fclose(f);
    table_fin(table);
    vtk_fin(vtk);
    csv_fin(csv);
    fprintf(stderr, "%s: %s\n", me, output);
  }
}
