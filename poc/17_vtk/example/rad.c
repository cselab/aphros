#include <stdio.h>
#include <stdlib.h>
#include <vtk.h>
#include <csv.h>
#include <table.h>
#include <math.h>

static double pi = 3.141592653589793;
#define Vcoef (4.0/3.0*pi)
static char me[] = "vtk/rad";

static double array_max(int, double*);
static double array_min(int, double*);

int
main(int argc, char **argv)
{
  int nt, nr, i, j, key, *flag;
  float *cl;
  double *cl_csv, *vf, *rad;
  struct VTK *vtk;
  struct CSV *csv;
  struct Table *table;
  char *path;
  FILE *f;
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
  for (i = 0; i < nt; i++)
      flag[i] = (cl[i] == -1 || cl[i] == 0);
  vtk_remove_tri(vtk, flag);
  free(flag);

  cl_csv = csv_field(csv, "cl");
  if (cl_csv == NULL) {
      fprintf(stderr, "%s: not field 'cl' in csv file\n", me);
      exit(2);
  }
  vf = csv_field(csv, "vf");
  if (cl_csv == NULL) {
      fprintf(stderr, "%s: not field 'vf' in csv file\n", me);
      exit(2);
  }
  nr = csv_nr(csv);
  fprintf(stderr, "min, max: %g %g\n", array_min(nr, vf),
	  array_max(nr, vf));
  table = table_ini(0, NULL, NULL);
  vtk_add(vtk, "rad", VTK_CELL, VTK_DOUBLE);
  rad = vtk_data(vtk, "rad");
  for (i = 0; i < nr; i++) {
      key = (int)cl_csv[i];
      table_put(table, &key, i);
  }

  nt = vtk_nt(vtk);
  for (i = 0; i < nt; i++)
      rad[i] = 0;

  for (i = 0; i < nt; i++) {
      key = (int)cl[i];
      if (key != 0) {
	  j = table_get(table, &key);
	  if (j != TABLE_EMPY)
	      rad[i] = pow(vf[j]/Vcoef, 1.0/3.0);
      }
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
