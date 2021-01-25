#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <csv.h>
#include <table.h>

enum { N = 999 };
static const char me[] = "split";

#include "util.h"

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
  fprintf(stderr, "%s -p prefix -f field [csv..]\n", me);
  exit(1);
}

#define GET(f, r)                                                  \
  if ((*r = csv_field(csv, f)) == NULL) {                          \
    fprintf(stderr, "%s: no field '%s' in '%s'\n", me, (f), name); \
    exit(2);                                                       \
  }

static const char* name;
struct Data {
  struct CSV* csv;
  struct Table* table;
  double* x;
  double* y;
  double* z;
  double* r;
  double* field;
};
static struct Data data_ini(const char*);
static int data_fin(struct Data*);
static int data_add(struct Data*, const char*, const int*);
static int dist(int, struct Data*, int, struct Data*, /**/ double*);

int main(int argc, char** argv) {
  char* Prefix;
  char output[N];
  double d;
  double dmin;
  FILE* file;
  int* array;
  int lmin;
  int m;
  int nb;
  int na;
  int nn;
  int i;
  int j;
  int l;
  int* new;
  int* prev;
  int* split;
  struct Data a;
  struct Data b;

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
  if (Prefix == NULL) {
    fprintf(stderr, "%s: prefix (-p) is not given\n", me);
    exit(1);
  }
  if (name == NULL) {
    fprintf(stderr, "%s: field (-f) is not given\n", me);
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
    for (i = 0; i < nb; i++) {
      split[i] = 0;
      prev[i] = -1;
    }
    /* nn: number of new */
    for (i = nn = 0; i < 2 * nb; i += 2)
      if (table_get(a.table, array[i], &j) == TABLE_EMPY)
        new[nn++] = array[i + 1];
    for (i = 0; i < nn; i++) {
      j = new[i];
      dmin = DBL_MAX;
      lmin = 0;
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
      if (table_get(b.table, (int)a.field[lmin], &m) != TABLE_EMPY) {
        split[m] = 1;
        split[j] = 2;
        prev[j] = prev[m] = (int)a.field[lmin];
        /*connect_add(j, &b);
           connect_add(m, &b); */
      } else {
        /*
           fprintf(stderr, "%s: prev disapeared: %s\n", me, *argv);
           fprintf(stderr, "%s: prev disapeared: lmin = %d, j = %d\n", me,
           lmin, j); */
      }
    }
    data_add(&b, "prev", prev);
    data_add(&b, "split", split);
    if (util_name(Prefix, *argv, output) != 0) {
      fprintf(stderr, "%s: util_name failed\n", me);
      exit(2);
    }
    if ((file = fopen(output, "w")) == NULL) {
      fprintf(stderr, "%s: fail to write to '%s'\n", me, output);
      exit(2);
    }
    if (csv_write(b.csv, file) != 0) {
      fprintf(stderr, "%s: csv_write failed\n", me);
      exit(2);
    }
    fclose(file);
    free(new);
    free(prev);
    free(split);
    free(array);
    data_fin(&a);
    a = b;
  }
  data_fin(&a);
}

static struct Data data_ini(const char* fname) {
  int nr, i;
  double* field;
  FILE* file;
  struct CSV* csv;
  struct Table* t;
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

static int data_fin(struct Data* q) {
  csv_fin(q->csv);
  return table_fin(q->table);
}

static int data_add(struct Data* q, const char* name, const int* a) {
  int n, i;
  double* f;
  struct CSV* csv;

  csv = q->csv;
  n = csv_nr(csv);
  if (csv_add(csv, name) != 0) {
    fprintf(stderr, "%s: csv_add failed\n", me);
    return 1;
  }
  if ((f = csv_field(csv, name)) == NULL) {
    fprintf(stderr, "%s: csv_field failed\n", me);
    return 1;
  }
  for (i = 0; i < n; i++)
    f[i] = a[i];
  return 0;
}

static int dist(int i, struct Data* a, int j, struct Data* b, /**/ double* p) {
  double d;
  double x;
  double y;
  double z;
  int na;
  int nb;

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
  // r = a->r[i] + b->r[j];
  d = sqrt(x * x + y * y + z * z);
  // d -= r;
  *p = d;
  return 0;
}
