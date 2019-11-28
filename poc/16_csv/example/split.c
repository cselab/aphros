#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <csv.h>
#include <table.h>

static const char me[] = "split";
#define	USED(x)		if(x);else{}
#define MALLOC(n, p)							\
    do {								\
	*(p) = malloc((n)*sizeof(**(p)));				\
	if (*(p) == NULL)  {						\
	    fprintf(stderr, "%s: alloc failed, n = %d", me, n);		\
	    exit(2);							\
	}								\
    } while(0)

static void
usg()
{
  fprintf(stderr, "%s -f field [csv..]\n", me);
  exit(1);
}

#define GET(f, r)							\
  if ((*r = csv_field(csv, f)) == NULL) {				\
    fprintf(stderr, "%s: no field '%s' in '%s'\n",			\
	    me, (f), name);						\
    exit(2);								\
  }

static const char *name;
struct Data {
  struct CSV *csv;
  struct Table *table;
  double *x;
  double *y;
  double *z;
  double *r;
  double *field;
};
static int data_fin(struct Data *);
static struct Data data_ini(const char *);
static int dist(int, struct Data *, int, struct Data *, /**/ double *);

int
main(int argc, char **argv)
{
  struct Data a;
  struct Data b;
  int *array;
  int *new;
  int *split;
  int *prev;
  int nb, na, nn, i, j, l;
  double d;
  double dmin;
  int lmin;
  int m;
  char *Prefix;

  USED(argc);
  name = Prefix = NULL;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
    case 'h':
      usg();
      break;
    case 'f':
      argv++;
      name = argv[0];
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
    array = table_array(b.table);
    if (array == NULL) {
      fprintf(stderr, "%s: table_array failed\n", me);
      exit(2);
    }
    nb = table_length(b.table);
    na = table_length(a.table);
    MALLOC(nb, &new);
    MALLOC(nb, &prev);
    MALLOC(nb, &split);
    /* nn: number of new */
    for (i = nn = 0; i < 2 * nb; i += 2)
      if (table_get(a.table, array[i], &j) == TABLE_EMPY)
        new[nn++] = array[i + 1];

    /* for (i = 0; i < k; i++) {
       j = new[i];
       printf("%s %d %g %g %g\n",
       *argv, (int) b.field[j], b.x[j], b.y[j], b.z[j]);
       } */
    for (i = 0; i < nn; i++) {
      j = new[i];
      dmin = DBL_MAX;
      for (l = 0; l < na; l++) {
        if (dist(j, &b, l, &a, &d) != 0) {
          fprintf(stderr, "%s: dist failed\n", me);
          exit(2);
        }
        if (d < dmin) {
          dmin = d;
          lmin = l;
        }
      }
      if (table_get(b.table, b.field[lmin], &m) != TABLE_EMPY) {
        split[m] = split[lmin] = 1;
        prev[m] = prev[lmin] = (int) b.field[lmin];
      } else {
        fprintf(stderr, "%s: prev disapeared: %s\n", me, *argv);
        fprintf(stderr, "%s: prev disapeared: lmin = %d, m = %d\n", me,
                lmin, m);
      }
    }
    free(new);
    free(prev);
    free(split);
    free(array);
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

  GET("x", &q.x);
  GET("y", &q.y);
  GET("z", &q.z);
  GET("r", &q.r);
  q.field = field;
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

static int
dist(int i, struct Data *a, int j, struct Data *b, /**/ double *p)
{
  double x;
  double y;
  double z;
  double r;
  double d;
  int na, nb;

  na = csv_nr(a->csv);
  nb = csv_nr(b->csv);
  if (i < 0 || i >= na) {
    fprintf(stderr, "%s: i=%d is not in [0, %d)\n", me, i, na);
    return 1;
  }
  if (j < 0 || j >= nb) {
    fprintf(stderr, "%s: j=%d is not in [0, %d)\n", me, j, nb);
    return 1;
  }
  x = a->x[i] - b->x[j];
  y = a->y[i] - b->y[j];
  z = a->z[i] - b->z[j];
  r = a->r[i] + b->r[j];
  d = sqrt(x * x + y * y + z * z);
  d -= r;
  *p = d;
  return 0;
}
