#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "err.h"
#include "memory.h"
#include "vtk.h"

char me[] = "vtk";
enum { N = 9999 };

#define FMT "%.20g"
#define LINE(s, f)				\
  do						\
    if (line(s, f) != 0) {			\
      MSG(("fail to read"));			\
      return NULL;				\
    }						\
  while(0)

#define SWAP(n, p) swap(n, sizeof(*(p)), p)
#define FILL(n, f, p)				\
    do {					\
	MALLOC(n, p);				\
	FREAD(n, (*(p)), f);			\
	SWAP(n, (*(p)));			\
    } while(0)
#define FREAD(n, p, f)							\
    do {								\
	if ((int)fread(p, sizeof(*(p)), (n), (f)) != (n))		\
	    ERR(("fread failed, n = %d", n));				\
    } while(0)
#define SIZE(a) (sizeof(a)/sizeof(*(a)))

static int swap(int, int, void *);
static int line(char *, FILE *);
static int line0(char *, FILE *);
static int eq(const char *, const char *);
static int rank2num(const char *, int *);
static int type2num(const char *, int *num, int *size);
static int location2num(const char *, int *);

static const char *LocationString[] = { "CELL_DATA", "POINT_DATA" };
static const int LocationEnum[] = { VTK_CELL, VTK_POINT };
static const char *TypeString[] = { "int", "float" };
static const int TypeEnum[] = { VTK_INT, VTK_FLOAT };
static const int TypeSize[] = { sizeof(int), sizeof(float) };
static const char *RankString[] = { "SCALARS", "VECTORS" };
static const int RankEnum[] = { VTK_SCALAR, VTK_VECTOR };

struct VTK *
vtk_read(FILE * f)
{
  struct VTK *q;
  double *x, *y, *z;
  void *p;
  float *r;
  int nv, nt, nf, n, m, i, j, nr, ifield, size;
  int *t, *t0, *t1, *t2;
  char s[N], rank[N], name[N], type[N], location[N];

  MALLOC(1, &q);
  LINE(s, f);
  if (!eq(s, "# vtk DataFile Version 2.0")) {
    MSG(("not a vtk file, got '%s'", s));
    return NULL;
  }
  LINE(s, f);                   // comments
  LINE(s, f);
  if (!eq(s, "BINARY")) {
    MSG(("expect 'BINARY', get '%s'", s));
    return NULL;
  }
  LINE(s, f);
  if (!eq(s, "DATASET POLYDATA")) {
    MSG(("expect 'DATASET POLYDATA', get '%s'", s));
    return NULL;
  }
  LINE(s, f);
  if (sscanf(s, "POINTS %d float", &nv) != 1) {
    MSG(("failt to parse: '%s'", s));
    return NULL;
  }
  FILL(3 * nv, f, &r);
  MALLOC(nv, &x);
  MALLOC(nv, &y);
  MALLOC(nv, &z);
  for (i = j = 0; i < nv; i++) {
    x[i] = r[j++];
    y[i] = r[j++];
    z[i] = r[j++];
  }
  FREE(r);
  LINE(s, f);
  if (sscanf(s, "POLYGONS %d %*d", &nt) != 1) {
    MSG(("fail to parse: '%s'", s));
    return NULL;
  }
  FILL(4 * nt, f, &t);
  MALLOC(nt, &t0);

  MALLOC(nt, &t1);
  MALLOC(nt, &t2);
  for (i = j = 0; i < nt; i++) {
    j++;
    t0[i] = t[j++];
    t1[i] = t[j++];
    t2[i] = t[j++];
  }
  FREE(t);
  nf = 0;
  if (line(s, f) != 0)
    goto end_data;
  for (;;) {
    sscanf(s, "%s", location);
    if (line(s, f) != 0)
      goto end_data;
    do {
      if (location2num(location, &q->location[nf]) != 0) {
        MSG(("unknown location '%s'", location));
        return NULL;
      }
      if (sscanf(s, "%s %s %s", rank, name, type) != 3)
        break;
      q->name[nf] = memory_strndup(name, N);
      if (eq(rank, "SCALARS")) {
        LINE(s, f);
        if (!eq(s, "LOOKUP_TABLE default")) {
          MSG(("expect 'LOOKUP_TABLE default', get '%s'", s));
          return NULL;
        }
      }
      if (rank2num(rank, &q->rank[nf]) != 0) {
        MSG(("unknown rank: '%s' for '%s'", rank, name));
        return NULL;
      }
      if (type2num(type, &q->type[nf], &size) != 0) {
        MSG(("unknown type: '%s' for '%s'", type, name));
        return NULL;
      }
      n = q->location[nf] == VTK_CELL ? nt : nv;
      m = n * q->rank[nf];
      MALLOC(m * size, &p);
      FREAD(m * size, p, f);
      swap(m, size, p);
      q->data[nf] = p;
      nf++;
    } while (line(s, f) == 0);
  }
end_data:
  q->nv = nv;
  q->nt = nt;
  q->nf = nf;
  q->x = x;
  q->y = y;
  q->z = z;
  q->t0 = t0;
  q->t1 = t1;
  q->t2 = t2;
  return q;
}

int
vtk_fin(struct VTK *q)
{
  int i, nf;

  nf = vtk_nf(q);
  for (i = 0; i < nf; i++) {
    FREE(q->name[i]);
    FREE(q->data[i]);
  }
  FREE(q->x);
  FREE(q->y);
  FREE(q->z);
  FREE(q->t0);
  FREE(q->t1);
  FREE(q->t2);
  FREE(q);
  return 0;
}

int
vtk_write(struct VTK *q, FILE * f)
{
  return 0;
}

static int
line(char *s, FILE * f)
{
  do
    if (line0(s, f) != 0)
      return 1;
  while (s[0] == '\0');
  return 0;
}

static int
line0(char *s, FILE * f)
{
  int n;

  if (fgets(s, N, f) == NULL)
    return 1;
  n = strlen(s);
  if (n > 0 && s[n - 1] == '\n')
    s[n - 1] = '\0';
  return 0;
}

int
vtk_nv(struct VTK *q)
{
  return q->nv;
}

int
vtk_nf(struct VTK *q)
{
  return q->nf;
}

int
vtk_nt(struct VTK *q)
{
  return q->nt;
}

void *
vtk_field(struct VTK *q, const char *name)
{
  int i, nf;

  nf = q->nf;
  for (i = 0; i < nf; i++)
    if (eq(name, q->name[i]))
      return q->data[i];
  return NULL;
}

static int
eq(const char *a, const char *b)
{
  return strncmp(a, b, N) == 0;
}

static int
swap(int n, int size, void *p0)
{
  int i;
  char *p, t;

  p = p0;
  while (n--) {
    for (i = 0; i < size / 2; i++) {
      t = p[i];
      p[i] = p[size - i - 1];
      p[size - i - 1] = t;
    }
    p += size;
  }
  return 0;
}

static int
location2num(const char *s, int *p)
{
  int i, n;

  n = SIZE(LocationString);
  for (i = 0; i < n; i++)
    if (eq(s, LocationString[i])) {
      *p = LocationEnum[i];
      return 0;
    }
  return 1;
}

static int
rank2num(const char *s, int *p)
{
  int i, n;

  n = SIZE(RankString);
  for (i = 0; i < n; i++)
    if (eq(s, RankString[i])) {
      *p = RankEnum[i];
      return 0;
    }
  return 1;
}

static int
type2num(const char *s, int *num, int *size)
{
  int i, n;

  n = SIZE(TypeString);
  for (i = 0; i < n; i++)
    if (eq(s, TypeString[i])) {
      *num = TypeEnum[i];
      *size = TypeSize[i];
      return 0;
    }
  return 1;
}
