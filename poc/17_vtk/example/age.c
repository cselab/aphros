#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <csv.h>
#include <table.h>

#define	USED(x)		if(x);else{}

static char me[] = "vtk/age";

int
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

int
main(int argc, char **argv)
{
  double a[9999], b[9999];
  struct Table *age;
  int na, nb, i, values, key, value, status;

  read(*++argv, &na, a);
  age = table_ini(100);
  for (i = 0; i < na; i++)
    table_put(age, a[i], 0);

  while (*++argv != NULL) {
    read(*argv, &nb, b);
    for (i = 0; i < nb; i++) {
      key = b[i];
      status = table_get(age, key, &value);
      if (status != TABLE_EMPY) {
        value += 1;
        value = -value;
        table_put(age, key, value);
      } else
        table_put(age, key, 0);
    }
    for (i = 0; i < nb; i++) {
      key = b[i];
      status = table_get(age, key, &value);
      if (status != TABLE_EMPY && value > 0)
        table_remove(age, key);
    }
    for (i = 0; i < nb; i++) {
      key = b[i];
      status = table_get(age, key, &value);
      assert(status != TABLE_EMPY);
      value = -value;
      table_put(age, key, value);
    }
  }

  for (i = 0; i < nb; i++) {
    key = b[i];
    status = table_get(age, key, &value);
    assert(status != TABLE_EMPY);
    printf("%d %d\n", key, value);
  }
  table_fin(age);
  USED(argc);
}
