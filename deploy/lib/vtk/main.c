#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "err.h"
#include "line.h"
#include "memory.h"
#include "vtk.h"

static char me[] = "vtk";
enum { N = 999 };

#define FMT "%.20g"
#define LINE(s, f)             \
  do                           \
    if (line_get(s, f) != 0) { \
      MSG(("fail to read"));   \
      goto fail;               \
    }                          \
  while (0)

#define SWAP(n, p) swap(n, sizeof(*(p)), p)
#define FILL(n, f, p)    \
  do {                   \
    MALLOC(n, p);        \
    FREAD(n, (*(p)), f); \
    SWAP(n, (*(p)));     \
  } while (0)
#define FREAD(n, p, f)                                \
  do {                                                \
    if ((int)fread(p, sizeof(*(p)), (n), (f)) != (n)) \
      ERR(("fread failed, n = %d", n));               \
  } while (0)
#define FWRITE(n, p, f)                                  \
  do {                                                   \
    if ((int)fwrite(p, sizeof(*(p)), (n), (f)) != (n)) { \
      MSG(("fail to write, n = %d", n));                 \
      return 1;                                          \
    }                                                    \
  } while (0)
#define SIZE(a) (sizeof(a) / sizeof(*(a)))

static const char* LocationString[] = {"CELL_DATA", "POINT_DATA"};
static const int LocationEnum[] = {VTK_CELL, VTK_POINT};
static const char* TypeString[] = {"int", "float", "double"};
static const int TypeEnum[] = {VTK_INT, VTK_FLOAT, VTK_DOUBLE};
static const int TypeSize[] = {sizeof(int), sizeof(float), sizeof(double)};
static const char* RankString[] = {"SCALARS", "VECTORS"};
static const int RankEnum[] = {VTK_SCALAR, VTK_VECTOR};

/* functions */
static int colormap(double, double, double, float*);
static double array_max(int, double*);
static double array_min(int, double*);
static int swap(int, int, void*);
static int eq(const char*, const char*);
static int rank2num(const char*, int*);
static int type2size(const char*, int* num, int* size);
static int location2num(const char*, int*);
static int num2location(int, const char**);
static int num2rank(int, const char**);
static int num2type(int, const char**);
static int num2size(int, int*);
static int remove0(int n, int size, void*, const int* a, int* new_size);

struct VTK* vtk_ini(
    int nv, const double* x, const double* y, const double* z, int nt,
    const int* t0, const int* t1, const int* t2) {
  int i;
  struct VTK* q;

  MALLOC(1, &q);
  MALLOC(nv, &q->x);
  MALLOC(nv, &q->y);
  MALLOC(nv, &q->z);
  MALLOC(nt, &q->t0);
  MALLOC(nt, &q->t1);
  MALLOC(nt, &q->t2);

  for (i = 0; i < nv; i++) {
    q->x[i] = x[i];
    q->y[i] = y[i];
    q->z[i] = z[i];
  }

  for (i = 0; i < nt; i++) {
    q->t0[i] = t0[i];
    q->t1[i] = t1[i];
    q->t2[i] = t2[i];
  }

  q->nv = nv;
  q->nt = nt;
  q->nf = 0;

  return q;
}

struct VTK* vtk_read(FILE* f) {
  char location[N];
  char name[N];
  char rank[N];
  char s[N];
  char type[N];
  double* x;
  double* y;
  double* z;
  float* r;
  int i;
  int j;
  int m;
  int n;
  int nf;
  int np0;
  int nt;
  int nv;
  int size;
  int* t;
  int* t0;
  int* t1;
  int* t2;
  struct VTK* q;
  void* p;

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
  LINE(s, f); // comments
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
  if (line_get(s, f) != 0) goto end_polygons;
  sscanf(s, "%[^\t ] %d float", name, &nv);
  if (!eq(name, "POINTS")) goto end_polygons;
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
  if (line_get(s, f) != 0) goto end_data;
  sscanf(s, "%[^\t ] %d %d", name, &nt, &np0);
  if (!eq(name, "POLYGONS")) goto end_polygons;
  if (4 * nt != np0) {
    MSG(("line: %s", s));
    MSG(("only triangles are supported"));
    goto fail;
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
  if (line_get(s, f) != 0) goto end_data;
end_polygons:
  for (;;) {
    sscanf(s, "%s", location);
    if (line_get(s, f) != 0) goto end_data;
    do {
      if (location2num(location, &q->location[nf]) != 0) {
        MSG(("unknown location '%s'", location));
        goto fail;
      }
      if (sscanf(s, "%[^\t ] %[^\t ] %[^\t ]", rank, name, type) != 3) break;
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
  line_fin();
  return q;
fail:
  line_fin();
  return NULL;
}

int vtk_fin(struct VTK* q) {
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

int vtk_write(struct VTK* q, FILE* f) {
  int nv, nt, nf, n, m, i, j, current, size, status;
  double *x, *y, *z;
  int *t0, *t1, *t2, *t;
  float* r;
  void* p;
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
    if (q->rank[i] == VTK_SCALAR) fprintf(f, "LOOKUP_TABLE default\n");
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

int vtk_off_color(struct VTK* q, const char* name, FILE* f) {
  int nt, nv, cnt[3];
  int *t0, *t1, *t2, i, j, index, ibuf[5], *di;
  double *x, *y, *z, *field, hi, lo, *dd;
  float *r, fbuf[3], *df;

  index = vtk_index(q, name);
  if (index == -1) {
    MSG(("no field '%s'", name));
    return 1;
  }
  if (q->rank[index] != 1) {
    MSG(("wrong rank=%d for '%s'", q->rank[index], name));
    return 1;
  }
  if (q->location[index] != VTK_CELL) {
    MSG(("wrong location for '%s'", name));
    return 1;
  }

  nt = vtk_nt(q);
  MALLOC(nt, &field);
  switch (q->type[index]) {
    case VTK_FLOAT:
      df = q->data[index];
      for (i = 0; i < nt; i++)
        field[i] = df[i];
      break;
    case VTK_DOUBLE:
      dd = q->data[index];
      for (i = 0; i < nt; i++)
        field[i] = dd[i];
      break;
    case VTK_INT:
      di = q->data[index];
      for (i = 0; i < nt; i++)
        field[i] = di[i];
      break;
    default:
      MSG(("wrong type for '%s'", name));
      return 1;
  }

  lo = array_min(nt, field);
  hi = array_max(nt, field);
  nv = vtk_nv(q);
  x = q->x;
  y = q->y;
  z = q->z;
  t0 = q->t0;
  t1 = q->t1;
  t2 = q->t2;
  fputs("OFF BINARY\n", f);

  cnt[0] = nv;
  cnt[1] = nt;
  cnt[2] = 0;
  SWAP(3, cnt);
  FWRITE(3, cnt, f);

  MALLOC(3 * nv, &r);
  for (i = j = 0; i < nv; i++) {
    r[j++] = x[i];
    r[j++] = y[i];
    r[j++] = z[i];
  }
  SWAP(3 * nv, r);
  FWRITE(3 * nv, r, f);
  FREE(r);

  for (i = 0; i < nt; i++) {
    j = 0;
    ibuf[j++] = 3;
    ibuf[j++] = t0[i];
    ibuf[j++] = t1[i];
    ibuf[j++] = t2[i];
    ibuf[j++] = 3;
    SWAP(j, ibuf);
    FWRITE(j, ibuf, f);
    colormap(field[i], lo, hi, fbuf);
    SWAP(3, fbuf);
    FWRITE(3, fbuf, f);
  }
  FREE(field);
  return 0;
}

int vtk_off(struct VTK* q, FILE* f) {
  int nt, nv, cnt[3];
  int *t, *t0, *t1, *t2, i, j;
  double *x, *y, *z;
  float* r;

  nt = vtk_nt(q);
  nv = vtk_nv(q);
  x = q->x;
  y = q->y;
  z = q->z;
  t0 = q->t0;
  t1 = q->t1;
  t2 = q->t2;
  fputs("OFF BINARY\n", f);

  cnt[0] = nv;
  cnt[1] = nt;
  cnt[2] = 0;
  SWAP(3, cnt);
  FWRITE(3, cnt, f);

  MALLOC(3 * nv, &r);
  for (i = j = 0; i < nv; i++) {
    r[j++] = x[i];
    r[j++] = y[i];
    r[j++] = z[i];
  }
  SWAP(3 * nv, r);
  FWRITE(3 * nv, r, f);
  FREE(r);

  MALLOC(5 * nt, &t);
  for (i = j = 0; i < nt; i++) {
    t[j++] = 3;
    t[j++] = t0[i];
    t[j++] = t1[i];
    t[j++] = t2[i];
    t[j++] = 0;
  }
  SWAP(5 * nt, t);
  FWRITE(5 * nt, t, f);
  FREE(t);

  return 0;
}

int vtk_remove_tri(struct VTK* q, const int* a) {
  int nt, nf, m, i, rank, size;
  void* data;

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

int vtk_remove_orphan(struct VTK* q) {
  int* a;
  int* b;
  int i;
  int m;
  int nf;
  int nt;
  int nv;
  int rank;
  int size;
  int* t0;
  int* t1;
  int* t2;
  void* data;

  nv = vtk_nv(q);
  nt = vtk_nt(q);
  nf = vtk_nf(q);
  t0 = q->t0;
  t1 = q->t1;
  t2 = q->t2;
  MALLOC(nv, &a);
  MALLOC(nv, &b);
  for (i = 0; i < nv; i++)
    a[i] = 1;
  for (i = 0; i < nt; i++) {
    a[t0[i]] = 0;
    a[t1[i]] = 0;
    a[t2[i]] = 0;
  }
  for (m = i = 0; i < nv; i++)
    if (a[i] == 0) b[i] = m++;
  for (i = 0; i < nt; i++) {
    t0[i] = b[t0[i]];
    t1[i] = b[t1[i]];
    t2[i] = b[t2[i]];
  }
  remove0(nv, sizeof(*q->x), q->x, a, &m);
  remove0(nv, sizeof(*q->y), q->y, a, &m);
  remove0(nv, sizeof(*q->z), q->z, a, &m);
  for (i = 0; i < nf; i++)
    if (q->location[i] == VTK_POINT) {
      rank = q->rank[i];
      num2size(q->type[i], &size);
      data = q->data[i];
      remove0(nv, size * rank, data, a, &m);
    }
  q->nv = m;
  FREE(a);
  FREE(b);
  return 0;
}

int vtk_nv(struct VTK* q) {
  return q->nv;
}

int vtk_nf(struct VTK* q) {
  return q->nf;
}

int vtk_nt(struct VTK* q) {
  return q->nt;
}

int vtk_index(struct VTK* q, const char* name) {
  int i, nf;

  nf = q->nf;
  for (i = 0; i < nf; i++)
    if (eq(name, q->name[i])) {
      return i;
    }
  return -1;
}

void* vtk_data(struct VTK* q, const char* name) {
  int i;

  i = vtk_index(q, name);
  if (i == -1)
    return NULL;
  else
    return q->data[i];
}

int vtk_remove(struct VTK* q, const char* name) {
  int i, j, nf;

  nf = vtk_nf(q);
  j = vtk_index(q, name);
  if (j == -1) return 1;
  FREE(q->name[j]);
  FREE(q->data[j]);
  for (i = j; i < nf - 1; i++) {
    if (i > j) FREE(q->name[i]);
    q->name[i] = memory_strndup(q->name[i + 1], N);
    q->location[i] = q->location[i + 1];
    q->type[i] = q->type[i + 1];
    q->rank[i] = q->rank[i + 1];
    q->data[i] = q->data[i + 1];
  }
  q->nf--;
  return 0;
}

int vtk_add(struct VTK* q, const char* name, int location, int type) {
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

static int eq(const char* a, const char* b) {
  return strncmp(a, b, N) == 0;
}

static int swap(int n, int size, void* p0) {
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

static int location2num(const char* s, int* p) {
  int i, n;

  n = SIZE(LocationString);
  for (i = 0; i < n; i++)
    if (eq(s, LocationString[i])) {
      *p = LocationEnum[i];
      return 0;
    }
  return 1;
}

static int rank2num(const char* s, int* p) {
  int i, n;

  n = SIZE(RankString);
  for (i = 0; i < n; i++)
    if (eq(s, RankString[i])) {
      *p = RankEnum[i];
      return 0;
    }
  return 1;
}

static int type2size(const char* s, int* num, int* size) {
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

static int num2type(int p, const char** s) {
  int i, n;

  n = SIZE(TypeString);
  for (i = 0; i < n; i++)
    if (p == TypeEnum[i]) {
      *s = TypeString[i];
      return 0;
    }
  return 1;
}

static int num2rank(int p, const char** s) {
  int i, n;

  n = SIZE(RankString);
  for (i = 0; i < n; i++)
    if (p == RankEnum[i]) {
      *s = RankString[i];
      return 0;
    }
  return 1;
}

static int num2location(int p, const char** s) {
  int i, n;

  n = SIZE(LocationString);
  for (i = 0; i < n; i++)
    if (p == LocationEnum[i]) {
      *s = LocationString[i];
      return 0;
    }
  return 1;
}

static int num2size(int p, int* s) {
  int i, n;

  n = SIZE(TypeEnum);
  for (i = 0; i < n; i++)
    if (p == TypeEnum[i]) {
      *s = TypeSize[i];
      return 0;
    }
  return 1;
}

static int remove0(int n, int size, void* pv, const int* a, int* pM) {
  int M, i, j, k, l;
  char* p;

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

static double array_min(int n, double* a) {
  int i;
  double m;

  if (n == 0) return 0;

  m = a[0];
  for (i = 1; i < n; i++)
    if (a[i] < m) m = a[i];
  return m;
}

static double array_max(int n, double* a) {
  int i;
  double m;

  if (n == 0) return 0;
  m = a[0];
  for (i = 1; i < n; i++)
    if (a[i] > m) m = a[i];
  return m;
}

static int colormap(double v, double l, double h, /**/ float p[3]) {
  float R, G, B;

  if (v < l) v = l;
  if (v > h) v = h;

  if (l != h)
    v = 4 * (v - l) / (h - l);
  else
    v = 0;

  R = 0;
  G = B = 1;
  if (v < 1)
    G = v;
  else if (v < 2)
    B = 2 - v;
  else if (v < 3) {
    R = v - 2;
    B = 0;
  } else {
    R = 1;
    G = 4 - v;
    B = 0;
  }

  p[0] = R;
  p[1] = G;
  p[2] = B;
  return 0;
}
