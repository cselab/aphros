#include <stdio.h>
#include <stdlib.h>
#include <csv.h>
#include <table.h>

static const char me[] = "split";

#define	USED(x)		if(x);else{}

static void
usg()
{
  fprintf(stderr, "%s -f field [csv..]\n", me);
  exit(1);
}

static const char *name;
struct Data {
  struct CSV *csv;
  struct Table *table;
  double *x;
  double *y;
  double *z;
  double *r;
};
static int data_fin(struct Data *);
static struct Data data_ini(const char *);

int
main(int argc, char **argv)
{
  struct Data a;
  struct Data b;
  int x, *y, n, i;

  USED(argc);
  name = NULL;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
    case 'h':
      usg();
      break;
    case 'f':
      argv++;
      name = argv[0];
      break;
    default:
      fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
      exit(1);
    }
  if (name == NULL) {
    fprintf(stderr, "%s: -f is not set\n", me);
    exit(1);
  }
  if (*argv == NULL) {
    fprintf(stderr, "%s: csv file is not given\n", me);
    exit(1);
  }

  a = data_ini(*argv);
  while (*++argv != NULL) {
      b = data_ini(*argv);
      y = table_array(b.table);
      n = table_length(b.table);
      for (i = 0; i < 2*n; i += 2) {
	  if (table_get(a.table, y[i], &x) == TABLE_EMPY)
	      printf("%s %d\n", *argv, y[i]);
      }
      free(y);
      data_fin(&a);
      a = b;
  }
  data_fin(&a);
}

static struct Data
data_ini(const char *fname)
{
  int nr, i;
  double *field;
  FILE *file;
  struct CSV *csv;
  struct Table *t;
  struct Data q;

  if ((file = fopen(fname, "r")) == NULL) {
    fprintf(stderr, "%s: fail to open '%s'\n", me, fname);
    exit(1);
  }
  if ((csv = csv_read(file)) == NULL) {
    fprintf(stderr, "%s: not a cvs file '%s'\n", me, fname);
    exit(1);
  }
  fclose(file);
  if ((field = csv_field(csv, name)) == NULL) {
    fprintf(stderr, "%s: no field '%s' in file '%s'\n", me, name, fname);
    exit(1);
  }
  nr = csv_nr(csv);
  t = table_ini(9999);
  for (i = 0; i < nr; i++) {
    table_put(t, field[i], i);
  }

  q.x = csv_field(csv, "x");
  q.y = csv_field(csv, "y");
  q.z = csv_field(csv, "z");
  q.r = csv_field(csv, "r");
  q.table = t;
  q.csv = csv;
  
  return q;
}

static int
data_fin(struct Data *q)
{
  csv_fin(q->csv);
  return table_fin(q->table);
}
