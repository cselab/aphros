#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <csv.h>
#include <table.h>

#define	USED(x)		if(x);else{}

static char me[] = "vtk/age";

static int read(const char*, int *, double*);
static int write(const char*, int, const double *, const double *);

int
main(int argc, char **argv)
{
  double cl[9999], age[9999];
  struct Table *tbl;
  int n, i, key, value, status, cnt;
  char name[999];

  cnt = 0;
  read(*++argv, &n, cl);
  tbl = table_ini(100);
  for (i = 0; i < n; i++) {
    age[i] = 0;
    table_put(tbl, cl[i], 0);
  }
  sprintf(name, "%05d.csv", cnt++);
  write(name, n, cl, age);

  while (*++argv != NULL) {
    read(*argv, &n, cl);
    for (i = 0; i < n; i++) {
      key = cl[i];
      status = table_get(tbl, key, &value);
      if (status != TABLE_EMPY) {
        value += 1;
        value = -value;
        table_put(tbl, key, value);
      } else
        table_put(tbl, key, 0);
    }
    for (i = 0; i < n; i++) {
      key = cl[i];
      status = table_get(tbl, key, &value);
      if (status != TABLE_EMPY && value > 0)
        table_remove(tbl, key);
    }
    for (i = 0; i < n; i++) {
      key = cl[i];
      status = table_get(tbl, key, &value);
      assert(status != TABLE_EMPY);
      value = -value;
      table_put(tbl, key, value);
      age[i] = value;
    }
    sprintf(name, "%05d.csv", cnt++);
    write(name, n, cl, age);
  }

  table_fin(tbl);
  USED(argc);
}

static int
read(const char *name, int *pn, double *a)
{
  int i, n;
  FILE *f;
  struct CSV *csv;
  double *field;

  if (name == NULL) {
    fprintf(stderr, "%s: name is empty\n", me);
    exit(2);
  }

  f = fopen(name, "r");
  if (f == NULL) {
    fprintf(stderr, "%s: fail to read '%s'\n", me, name);
    exit(2);
  }
  csv = csv_read(f);
  fclose(f);
  if (csv == NULL) {
    fprintf(stderr, "%s: fail to parse '%s'\n", me, name);
    exit(2);
  }
  field = csv_field(csv, "cl");
  if (field == NULL) {
    fprintf(stderr, "%s: no field '%s'\n", me, "cl");
    exit(2);
  }
  n = csv_nr(csv);
  for (i = 0; i < n; i++)
    a[i] = field[i];
  *pn = n;
  csv_fin(csv);
  return 0;
}

static int
write(const char *name, int n, const double *cl, const double *age)
{
  int i;
  struct CSV *csv;
  FILE *f;
  double *cl0, *age0;
  
  csv = csv_ini(n);
  if (csv == NULL) {
    fprintf(stderr, "%s: allocate\n", me);
    exit(2);
  }
  
  csv_add(csv, "cl");
  csv_add(csv, "age");

  cl0 = csv_field(csv, "cl");
  age0 = csv_field(csv, "age");

  for (i = 0; i < n; i++) {
    cl0[i] = cl[i];
    age0[i] = age[i];
  }

  f = fopen(name, "w");
  if (f == NULL) {
    fprintf(stderr, "%s: fail to open '%s'\n", me, name);
    exit(2);
  }
  csv_write(csv, f);
  fclose(f);
  csv_fin(csv);
  return 0;
}
