#include <stdio.h>
#include <string.h>
#include "err.h"
#include "memory.h"
#include "vtk.h"

char me[] = "vtk";
enum { N = 9999 };

#define FMT "%.20g"
#define LINE(s, f)				\
    if (line(s, f) != 0)			\
	ERR(("fail to read"))
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

static int swap(int, int, void *);
static int line(char *, FILE *);
static int eq(const char *, const char *);

struct VTK *
vtk_read(FILE * f)
{
  struct VTK *q;
  char s[N];
  const char *field;
  double *x, *y, *z;
  float *r;
  int nv, nt, nf, i, j, nr, ifield, M;
  int *t, *t0, *t1, *t2;

  MALLOC(1, &q);
  LINE(s, f);
  if (!eq(s, "# vtk DataFile Version 2.0")) {
    MSG(("not a vtk file, got '%s'", s));
    return NULL;
  }
  LINE(s, f);                   // comments
  do
    LINE(s, f);
  while (s[0] == '\0');
  if (!eq(s, "BINARY")) {
    MSG(("expect 'BINARY', get '%s'", s));
    return NULL;
  }
  do
    LINE(s, f);
  while (s[0] == '\0');
  if (!eq(s, "DATASET POLYDATA")) {
    MSG(("expect 'DATASET POLYDATA', get '%s'", s));
    return NULL;
  }
  do
    LINE(s, f);
  while (s[0] == '\0');
  if (sscanf(s, "POINTS %d float", &nv) != 1) {
    MSG(("failt to parse: '%s'", s));
    return NULL;
  }
  MSG(("nv: %d", nv));
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
  do
      LINE(s, f);
  while (s[0] == '\0');
  if (sscanf(s, "POLYGONS %d %*d", &nt) != 1) {
	MSG(("fail to parse: '%s'", s));
	return NULL;
  }
  FILL(4 * nt, stdin, &t);
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
  
  M = 2;                        /* initial size */
  nf = 0;
  field = strtok(s, ",");
  while (field != NULL) {
    if (nf == VTK_MAX_NF)
      ERR(("nf == VTK_MAX_NF=%d", nf, VTK_MAX_NF));
    q->name[nf] = memory_strndup(field, N);
    MALLOC(M, &q->data[nf]);
    nf++;
    field = strtok(NULL, ",");
  }

  q->nv = nv;
  q->x = x;
  q->y = y;
  q->z = z;
  q->nt = nt;
  q->t0 = t0;
  q->t1 = t1;
  q->t2 = t2;
  return q;
}

int
vtk_fin(struct VTK *q)
{
  int i, nf;

  nf = q->nf;
  for (i = 0; i < nf; i++) {
    FREE(q->name[i]);
    FREE(q->data[i]);
  }
  FREE(q->x);
  FREE(q->y);
  FREE(q->z);
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
vtk_nt(struct VTK *q)
{
  return q->nt;
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
