#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <vtk.h>
#include <csv.h>
#include <table.h>

#define	USED(x)		if(x);else{}
static double pi = 3.141592653589793;

#define Vcoef (4.0/3.0*pi)
static char me[] = "vtk/rad";

static double array_max(int, double *);
static double array_min(int, double *);
static int array_max_arg(int, double *);
static int array_min_arg(int, double *);

int
main(int argc, char **argv)
{
  int status, nt, nr, i, j, key, *flag;
  float *cl;
  double *cl_csv, *vf, *rad;
  struct VTK *vtk;
  struct CSV *csv;
  struct Table *table;
  char *path;
  FILE *f;

  USED(argc);

  argv++;
  path = argv++[0];
  if (path == NULL) {
    fprintf(stderr, "%s: needs vtk file\n", me);
    exit(2);
  }
  f = fopen(path, "r");
  if (f == NULL) {
    fprintf(stderr, "%s: cannot open '%s'\n", me, path);
    exit(2);
  }
  vtk = vtk_read(f);
  if (vtk == NULL) {
    fprintf(stderr, "%s: fail to parse csv file '%s'\n", me, path);
    exit(2);
  }
  fclose(f);

  path = argv++[0];
  if (path == NULL) {
    fprintf(stderr, "%s: needs csv file\n", me);
    exit(2);
  }
  f = fopen(path, "r");
  if (f == NULL) {
    fprintf(stderr, "%s: cannot open '%s'\n", me, path);
    exit(2);
  }
  csv = csv_read(f);
  if (csv == NULL) {
    fprintf(stderr, "%s: fail to parse csv file '%s'\n", me, path);
    exit(2);
  }
  fclose(f);
  nt = vtk_nt(vtk);
  flag = malloc(nt * sizeof(*flag));
  cl = vtk_data(vtk, "cl");
  if (cl == NULL) {
    fprintf(stderr, "%s: not field 'cl' in vtk file\n", me);
    exit(2);
  }
  cl_csv = csv_field(csv, "cl");
  if (cl_csv == NULL) {
    fprintf(stderr, "%s: not field 'cl' in csv file\n", me);
    exit(2);
  }
  vf = csv_field(csv, "vf");
  if (vf == NULL) {
    fprintf(stderr, "%s: not field 'vf' in csv file\n", me);
    exit(2);
  }
  nr = csv_nr(csv);
  j = array_max_arg(nr, vf);
  for (i = 0; i < nt; i++) {
    key = (int) cl[i];
    flag[i] = (key == -1 || key == 0 || key == cl_csv[j]);
  }
  vtk_remove_tri(vtk, flag);
  free(flag);
  vtk_remove(vtk, "nn");
  vtk_remove(vtk, "c");
  vtk_remove(vtk, "l");

  table = table_ini(100);
  for (i = 0; i < nr; i++) {
    key = (int) cl_csv[i];
    table_put(table, key, i);
  }
  vtk_add(vtk, "rad", VTK_CELL, VTK_DOUBLE);
  rad = vtk_data(vtk, "rad");
  nt = vtk_nt(vtk);
  for (i = 0; i < nt; i++) {
    key = (int) cl[i];
    status = table_get(table, key, &j);
    if (status != TABLE_EMPY)
      rad[i] = pow(vf[j] / Vcoef, 1.0 / 3.0);
    else
      rad[i] = 0;
  }
  vtk_write(vtk, stdout);
  table_fin(table);
  vtk_fin(vtk);
  csv_fin(csv);
}

static double
array_min(int n, double *a)
{
  int i;
  double m;

  m = a[0];
  for (i = 1; i < n; i++)
    if (a[i] < m)
      m = a[i];
  return m;
}

static double
array_max(int n, double *a)
{
  int i;
  double m;

  m = a[0];
  for (i = 1; i < n; i++)
    if (a[i] > m)
      m = a[i];
  return m;
}

static int
array_min_arg(int n, double *a)
{
  int i, j;
  double m;

  j = 0;
  m = a[j];
  for (i = 1; i < n; i++)
    if (a[i] < m) {
      j = i;
      m = a[j];
    }
  return j;
}

static int
array_max_arg(int n, double *a)
{
  int i, j;
  double m;

  j = 0;
  m = a[j];
  for (i = 1; i < n; i++)
    if (a[i] > m) {
      j = i;
      m = a[j];
    }
  return j;
}
