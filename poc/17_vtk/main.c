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
	ERR(("fail to read"));
static int line(char *, FILE *);
static int eq(const char *, const char *);

struct VTK *
vtk_read(FILE * f)
{
  struct VTK *q;
  char s[N];
  const char *field;
  int i, nf, nr, ifield, M;

  MALLOC(1, &q);
  LINE(s, f);
  if (!eq(s, "# vtk DataFile Version 2.0")) {
      MSG(("not a vtk file, got '%s'", s));
    return NULL;
  }
  
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

  nr = 0;
  while (line(s, f) == 0) {
    ifield = 0;
    field = strtok(s, ",");
    while (field != NULL) {
      if (ifield == nf)
        ERR(("ifield == nf=%d, nr=%d", nf, nr));
      if (ifield == 0 && nr == M) {
        M *= 2;
        for (i = 0; i < nf; i++)
          REALLOC(M, &q->data[i]);
      }
      ifield++;
      field = strtok(NULL, ",");
    }
    if (ifield != nf)
      ERR(("ifield=%d != nf=%d, nr=%d", ifield, nf, nr));
    nr++;
  }

  q->nr = nr;
  q->nf = nf;
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

static int
eq(const char *a, const char *b)
{
  return strncmp(a, b, N) == 0;
}
