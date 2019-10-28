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
static int type2size(const char *, int *num, int *size);
static int location2num(const char *, int *);

static int num2location(int, const char **);
static int num2rank(int, const char **);
static int num2type(int, const char **);
static int num2size(int, int *);
static int remove0(int n, int size, void *, const int *a, int *new_size);

struct VTK *
vtk_read(FILE * f)
{
  struct VTK *q;
  double *x, *y, *z;
  void *p;
  float *r;
  int nv, nt, nf, n, m, i, j, size;
  int *t, *t0, *t1, *t2;
  char s[N], rank[N], name[N], type[N], location[N];

  nv = nt = nf = 0;
  MALLOC(1, &q);
  line_ini();
  x = y = z = NULL;
  t0 = t1 = t2 = NULL;
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
      if (type2size(type, &q->type[nf], &size) != 0) {
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
  int nv, nt, nf, n, m, i, j, current, size, status;
  double *x, *y, *z;
  int *t0, *t1, *t2, *t;
  float *r;
  void *p;
  const char *rank, *type, *location;

  nv = q->nv;
  nt = q->nt;
  x = q->x;
  y = q->y;
  z = q->z;
  t0 = q->t0;
  t1 = q->t1;
  t2 = q->t2;
  if (fputs("# vtk DataFile Version 2.0\n", f) == EOF) {
    MSG(("fail to write"));
    return 1;
  }
  fprintf(f, "%s\n", me);
  fputs("BINARY\n", f);
  fputs("DATASET POLYDATA\n", f);
  fprintf(f, "POINTS %d float\n", nv);
  if (nv > 0) {
    MALLOC(3 * nv, &r);
    for (i = j = 0; i < nv; i++) {
      r[j++] = x[i];
      r[j++] = y[i];
      r[j++] = z[i];
    }
    SWAP(3 * nv, r);
    FWRITE(3 * nv, r, f);
    FREE(r);
  }
  fprintf(f, "POLYGONS %d %d\n", nt, 4 * nt);
  if (nt > 0) {
    MALLOC(4 * nt, &t);
    for (i = j = 0; i < nt; i++) {
      t[j++] = 3;
      t[j++] = t0[i];
      t[j++] = t1[i];
      t[j++] = t2[i];
    }
    SWAP(4 * nt, t);
    FWRITE(4 * nt, t, f);
    FREE(t);
  }

  size = current = n = -1;
  type = rank = NULL;
  nf = vtk_nf(q);
  for (i = 0; i < nf; i++) {
    if (q->location[i] != current) {
      current = q->location[i];
      if (current == VTK_CELL)
        n = nt;
      else if (current == VTK_POINT)
        n = nv;
      else {
        MSG(("unknown location: %d", current));
        goto end;
      }
      num2location(q->location[i], &location);
      fprintf(f, "%s %d\n", location, n);
    }
    num2rank(q->rank[i], &rank);
    num2type(q->type[i], &type);
    status = num2size(q->type[i], &size);
    if (status != 0) {
      MSG(("unknown type: %d", type));
      goto end;
    }
    fprintf(f, "%s %s %s\n", rank, q->name[i], type);
    if (q->rank[i] == VTK_SCALAR)
      fprintf(f, "LOOKUP_TABLE default\n");
    m = n * q->rank[i];
    MALLOC(m * size, &p);
    memory_memcpy(p, q->data[i], m * size);
    swap(m, size, p);
    FWRITE(m * size, p, f);
    FREE(p);
  }
end:
  return 0;
}

int
vtk_remove_tri(struct VTK *q, const int *a)
{
  int nt, nf, m, i, rank, size;
  void *data;

  nt = vtk_nt(q);
  nf = vtk_nf(q);

  size = sizeof(*q->t0);
  remove0(nt, size, q->t0, a, &m);
  remove0(nt, size, q->t1, a, &m);
  remove0(nt, size, q->t2, a, &m);
  for (i = 0; i < nf; i++)
    if (q->location[i] == VTK_CELL) {
      rank = q->rank[i];
      num2size(q->type[i], &size);
      data = q->data[i];
      remove0(nt, size * rank, data, a, &m);
    }
  q->nt = m;
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

int
vtk_index(struct VTK *q, const char *name, int *p)
{
  int i, nf;

  nf = q->nf;
  for (i = 0; i < nf; i++)
      if (eq(name, q->name[i])) {
	*p = i;
	return 0;
      }
  return 1;
}

void *
vtk_data(struct VTK *q, const char *name)
{
  int i, status;
  status = vtk_index(q, name, &i);
  if (status != 0)
      return NULL;
  else
      return q->data[i];
}

int
vtk_add(struct VTK * q, const char *name, int location, int type)
{
    int nt, nv, nf, n, size, rank;

    nf = vtk_nf(q);
    nv = vtk_nv(q);
    nt = vtk_nt(q);
    if (nf == VTK_MAX_NF) {
	MSG(("nf=%d == VTK_MAX_NF", nf));
	return 1;
    }

    rank = VTK_SCALAR;
    q->location[nf] = location;
    q->type[nf] = type;
    q->rank[nf] = rank;
    q->name[nf] = memory_strndup(name, N);
    num2size(type, &size);
    if (location == VTK_CELL)
	n = nt;
    else if (location == VTK_POINT)
	n = nv;
    else
	MSG(("unknown location: %d", location));
    MALLOC(n * size * rank, &q->data[nf]);
    nf++;
    q->nf = nf;
    return 0;
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
type2size(const char *s, int *num, int *size)
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

static int
num2type(int p, const char **s)
{
  int i, n;

  n = SIZE(TypeString);
  for (i = 0; i < n; i++)
    if (p == TypeEnum[i]) {
      *s = TypeString[i];
      return 0;
    }
  return 1;
}

static int
num2rank(int p, const char **s)
{
  int i, n;

  n = SIZE(RankString);
  for (i = 0; i < n; i++)
    if (p == RankEnum[i]) {
      *s = RankString[i];
      return 0;
    }
  return 1;
}

static int
num2location(int p, const char **s)
{
  int i, n;

  n = SIZE(LocationString);
  for (i = 0; i < n; i++)
    if (p == LocationEnum[i]) {
      *s = LocationString[i];
      return 0;
    }
  return 1;
}

static int
num2size(int p, int *s)
{
  int i, n;

  n = SIZE(TypeEnum);
  for (i = 0; i < n; i++)
    if (p == TypeEnum[i]) {
      *s = TypeSize[i];
      return 0;
    }
  return 1;
}

static int
remove0(int n, int size, void *pv, const int *a, int *pM)
{
  int M, i, j, k, l;
  char *p;

  p = pv;
  M = 0;
  for (i = k = l = 0; i < n; i++)
    if (a[i] == 0) {
      for (j = 0; j < size; j++)
        p[k++] = p[size * i + j];
      M++;
    }
  *pM = M;
  return 0;
}
