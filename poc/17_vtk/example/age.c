#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include <csv.h>
#include <table.h>

#define	USED(x)		if(x);else{}

enum { N = 999, M = 99999 };
static char me[] = "vtk/age";

static int read(const char *, int *, double *);
static int write(const char *, int, const double *, const double *);
static int digits(const char *, char *);

static void
usg(void)
{
  fprintf(stderr, "%s [-o pattern] [file.csv ..]\n", me);
  exit(1);
}

int
main(int argc, char **argv)
{
  double *cl, *age;
  struct Table *tbl;
  int n, i, key, value, status;
  char dig[N], name[N];
  char *pattern;

  USED(argc);
  pattern = NULL;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
    case 'h':
      usg();
      break;
    case 'p':
      if ((pattern = *argv++) == NULL) {
        fprintf(stderr, "%s: -p needs an argument\n", me);
        exit(2);
      }
      break;
    default:
      fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
      exit(1);
    }
  if (*argv == NULL) {
    fprintf(stderr, "%s: no files given\n", me);
    exit(1);
  }
  if (pattern == NULL) {
    fprintf(stderr, "%s: pattern (-p) is not given\n", me);
    exit(1);
  }
  cl = malloc(M * sizeof(*cl));
  if (cl == NULL) {
    fprintf(stderr, "alloc failed\n");
    exit(2);
  }
  age = malloc(M * sizeof(*age));
  if (age == NULL) {
    fprintf(stderr, "alloc failed\n");
    exit(2);
  }
  read(*argv, &n, cl);
  if (n > M) {
    fprintf(stderr, "n=%d > M=%d\n", n, M);
    exit(2);
  }
  tbl = table_ini(100);
  for (i = 0; i < n; i++) {
    age[i] = 0;
    table_put(tbl, cl[i], 0);
  }
  digits(*argv, dig);
  if (snprintf(name, N, "%s.csv", dig) < 0) {
    fprintf(stderr, "snprintf failed\n");
    exit(2);
  }
  write(name, n, cl, age);

  while (*++argv != NULL) {
    read(*argv, &n, cl);
    if (n > M) {
      fprintf(stderr, "n=%d > M=%d\n", n, M);
      exit(2);
    }
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
    digits(*argv, dig);
    if (snprintf(name, N, "%s.csv", dig) < 0) {
      fprintf(stderr, "snprintf failed\n");
      exit(2);
    }
    write(name, n, cl, age);
  }

  table_fin(tbl);
  free(cl);
  free(age);
  return 0;
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

static int
digits(const char *s, char *ans)
{
  int i, j;

  for (i = j = 0; s[i] != '\0'; i++)
    if (s[i] == '/')
      j = i + 1;
  for (; s[j] != '\0'; j++)
    if (isdigit(s[j]))
      break;
  for (i = j; s[i] != '\0'; i++)
    if (!isdigit(s[i]))
      break;
    else
      *ans++ = s[i];
  *ans = '\0';
  return 0;
}