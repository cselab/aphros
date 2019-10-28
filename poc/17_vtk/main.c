#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "err.h"
#include "memory.h"
#include "line.h"
#include "vtk.h"

static char me[] = "vtk";
enum { N = 999 };

#define FMT "%.20g"
#define LINE(s, f)				\
  do						\
    if (line_get(s, f) != 0) {			\
      MSG(("fail to read"));			\
      goto fail;				\
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
#define FWRITE(n, p, f)							\
  do {									\
    if ((int)fwrite(p, sizeof(*(p)), (n), (f)) != (n)) {		\
      MSG(("failt to write, n = %d\n", n));				\
      return 1;								\
    }									\
  } while(0)
#define SIZE(a) (sizeof(a)/sizeof(*(a)))

static const char *LocationString[] = { "CELL_DATA", "POINT_DATA" };
static const int LocationEnum[] = { VTK_CELL, VTK_POINT };
static const char *TypeString[] = { "int", "float", "double" };
static const int TypeEnum[] = { VTK_INT, VTK_FLOAT, VTK_DOUBLE };
static const int TypeSize[] =
    { sizeof(int), sizeof(float), sizeof(double) };
static const char *RankString[] = { "SCALARS", "VECTORS" };
static const int RankEnum[] = { VTK_SCALAR, VTK_VECTOR };

/* functions */
static int swap(int, int, void *);
static int eq(const char *, const char *);
static int rank2num(const char *, int *);
static int type2num(const char *, int *num, int *size);
static int location2num(const char *, int *);

struct VTK *
vtk_read(FILE * f)
{
  struct VTK *q;
  double *x, *y, *z;
  void *p;
  float *r;
  int nv, nt, nf, n, m, i, j, nr, ifield, size;
  int *t, *t0, *t1, *t2;
  char s[N], u[N], rank[N], name[N], type[N], location[N];

  nv = nt = nf = 0;
  MALLOC(1, &q);
  line_ini();

  LINE(s, f);
  if (!eq(s, "# vtk DataFile Version 2.0")) {
    MSG(("not a vtk file, got '%s'", s));
    goto fail;
  }
  LINE(s, f);                   // comments
  LINE(s, f);
  if (!eq(s, "BINARY")) {
    MSG(("expect 'BINARY', get '%s'", s));
    goto fail;
  }
  LINE(s, f);
  if (!eq(s, "DATASET POLYDATA")) {
    MSG(("expect 'DATASET POLYDATA', get '%s'", s));
    goto fail;
  }
  if (line_get(s, f) != 0)
    goto end_polygons;
  sscanf(s, "%s %d float", name, &nv);
  if (!eq(name, "POINTS"))
    goto end_polygons;
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
  if (line_get(s, f) != 0)
    goto end_data;
  sscanf(s, "%s %d %*d", name, &nt);
  if (!eq(name, "POLYGONS")) {
    MSG(("preved: %s", name));
    goto end_polygons;
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
  if (line_get(s, f) != 0)
    goto end_data;
end_polygons:
  for (;;) {
    sscanf(s, "%s", location);
    if (line_get(s, f) != 0)
      goto end_data;
    do {
      if (location2num(location, &q->location[nf]) != 0) {
        MSG(("unknown location '%s'", location));
        goto fail;
      }
      if (sscanf(s, "%s %s %s", rank, name, type) != 3)
        break;
      q->name[nf] = memory_strndup(name, N);
      if (eq(rank, "SCALARS")) {
        LINE(s, f);
        if (!eq(s, "LOOKUP_TABLE default")) {
          MSG(("expect 'LOOKUP_TABLE default', get '%s'", s));
          goto fail;
        }
      }
      if (rank2num(rank, &q->rank[nf]) != 0) {
        MSG(("unknown rank: '%s' for '%s'", rank, name));
        goto fail;
      }
      if (type2num(type, &q->type[nf], &size) != 0) {
        MSG(("unknown type: '%s' for '%s'", type, name));
        goto fail;
      }
      n = q->location[nf] == VTK_CELL ? nt : nv;
      m = n * q->rank[nf];
      MALLOC(m * size, &p);
      FREAD(m * size, p, f);
      swap(m, size, p);
      q->data[nf] = p;
      nf++;
    } while (line_get(s, f) == 0);
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
  //line_write(stderr);
  line_fin();
  return q;
fail:
  line_fin();
  return NULL;
}

int
vtk_fin(struct VTK *q)
{
  int i, nv, nt, nf;

  nv = vtk_nv(q);
  nt = vtk_nt(q);
  nf = vtk_nf(q);
  for (i = 0; i < nf; i++) {
    FREE(q->name[i]);
    FREE(q->data[i]);
  }
  if (nv > 0) {
    FREE(q->x);
    FREE(q->y);
    FREE(q->z);
  }
  if (nt > 0) {
    FREE(q->t0);
    FREE(q->t1);
    FREE(q->t2);
  }
  FREE(q);
  return 0;
}

int
vtk_write(struct VTK *q, FILE * f)
{
  int nv, nt, i, j;
  double *x, *y, *z;
  float *r;

  nv = q->nv;
  nt = q->nt;
  x = q->x;
  y = q->y;
  z = q->z;
  if (fputs("# vtk DataFile Version 2.0\n", f) == EOF) {
    MSG(("fail to write"));
    return 1;
  }
  fprintf(f, "%s\n", me);
  fputs("BINARY\n", f);
  fputs("DATASET POLYDATA\n", f);
  if (nv > 0)
    fprintf(f, "POINTS %d float\n", nv);
  MALLOC(3 * nv, &r);
  for (i = j = 0; i < nv; i++) {
    r[j++] = x[i];
    r[j++] = y[i];
    r[j++] = z[i];
  }
  SWAP(3 * nv, r);
  FWRITE(3 * nv, r, f);
  FREE(r);
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
