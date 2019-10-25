#include <stdio.h>
#include <string.h>
#include "err.h"
#include "memory.h"
#include "csv.h"

char me[] = "csv";
enum { N = 9999 };

#define LINE(s, f)				\
    if (line(s, f) != 0)			\
	ERR(("fail to read"));
static int line(char *, FILE *);
static int eq(const char *, const char *);

struct CSV *
csv_read(FILE * f)
{
  struct CSV *q;
  char s[N];
  const char *field;
  int i, nf, nr, ifield, M;

  MALLOC(1, &q);
  LINE(s, f);
  M = 2;                        /* initial size */
  nf = 0;
  field = strtok(s, ",");
  while (field != NULL) {
    if (nf == CSV_MAX_NF)
      ERR(("nf == CSV_MAX_NF=%d", nf, CSV_MAX_NF));
    MSG(("%s", field));
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
        ERR(("ifield == nf=%d, nr = %d", nf, nr));
      if (ifield == 0 && nr == M) {
        M *= 2;
        for (i = 0; i < nf; i++)
          REALLOC(M, &q->data[i]);
      }
      if (sscanf(field, "%lf", &q->data[ifield][nr]) != 1)
        ERR(("not a number '%s'", field));
      ifield++;
      field = strtok(NULL, ",");
    }
    nr++;
  }

  q->nr = nr;
  q->nf = nf;
  return q;
}

int
csv_fin(struct CSV *q)
{
  int i, nf;

  nf = q->nf;
  for (i = 0; i < nf; i++) {
    FREE(q->name[i]);
    FREE(q->data[i]);
  }
  FREE(q);
}

double *
csv_field(struct CSV *q, const char *name)
{
  int i, nf;

  nf = q->nf;
  for (i = 0; i < nf; i++)
    if (eq(name, q->name[i]))
      return q->data[i];
  return NULL;
}

int
csv_nf(struct CSV *q)
{
  return q->nf;
}

int
csv_nr(struct CSV *q)
{
  return q->nr;
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
