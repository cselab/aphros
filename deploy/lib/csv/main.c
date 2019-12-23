#include <stdio.h>
#include <string.h>

#include "csv.h"
#include "err.h"
#include "memory.h"

static char me[] = "csv";
enum { N = 9999 };

#define FMT "%.20g"
#define LINE(s, f) \
  if (line(s, f) != 0) ERR(("fail to read"));
static int line(char*, FILE*);
static int eq(const char*, const char*);

struct CSV* csv_ini(int nr) {
  struct CSV* q;
  MALLOC(1, &q);

  q->nr = nr;
  q->nf = 0;

  return q;
}

struct CSV* csv_read(FILE* f) {
  struct CSV* q;
  char s[N];
  const char* field;
  int i, nf, nr, ifield, M;

  MALLOC(1, &q);
  LINE(s, f);
  M = 2; /* initial size */
  nf = 0;
  field = strtok(s, ",");
  while (field != NULL) {
    if (nf == CSV_MAX_NF) {
      MSG(("nf == CSV_MAX_NF=%d", nf, CSV_MAX_NF));
      return NULL;
    }
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
      if (ifield == nf) {
        MSG(("ifield == nf=%d, nr=%d", nf, nr));
        return NULL;
      }
      if (ifield == 0 && nr == M) {
        M *= 2;
        for (i = 0; i < nf; i++)
          REALLOC(M, &q->data[i]);
      }
      if (sscanf(field, "%lf", &q->data[ifield][nr]) != 1) {
        MSG(("not a number '%s'", field));
        return NULL;
      }
      ifield++;
      field = strtok(NULL, ",");
    }
    if (ifield != nf) {
      MSG(("ifield=%d != nf=%d, nr=%d", ifield, nf, nr));
      return NULL;
    }
    nr++;
  }

  q->nr = nr;
  q->nf = nf;
  return q;
}

int csv_fin(struct CSV* q) {
  int i, nf;

  nf = q->nf;
  for (i = 0; i < nf; i++) {
    FREE(q->name[i]);
    FREE(q->data[i]);
  }
  FREE(q);
  return 0;
}

double* csv_field(struct CSV* q, const char* name) {
  int i, nf;

  nf = q->nf;
  for (i = 0; i < nf; i++)
    if (eq(name, q->name[i])) return q->data[i];
  return NULL;
}

int csv_nf(struct CSV* q) {
  return q->nf;
}

int csv_nr(struct CSV* q) {
  return q->nr;
}

int csv_add(struct CSV* q, const char* name) {
  int nr, nf;
  double* data;

  nr = csv_nr(q);
  nf = csv_nf(q);
  MALLOC(nr, &data);
  q->data[nf] = data;
  q->name[nf] = memory_strndup(name, N);
  nf++;
  q->nf = nf;
  return 0;
}

static int write(struct CSV* q, char sep, FILE* f) {
  int nf, nr, i, j;

  nf = csv_nf(q);
  nr = csv_nr(q);
  for (i = 0; i < nf; i++) {
    if (i > 0) fputc(sep, f);
    if (fputs(q->name[i], f) == EOF) {
      MSG(("fail to write"));
      return 1;
    }
  }
  fputc('\n', f);
  for (j = 0; j < nr; j++) {
    for (i = 0; i < nf; i++) {
      if (i > 0) fputc(sep, f);
      fprintf(f, FMT, q->data[i][j]);
    }
    fputc('\n', f);
  }
  return 0;
}

int csv_write(struct CSV* q, FILE* f) {
  return write(q, ',', f);
}

int csv_write_space(struct CSV* q, FILE* f) {
  return write(q, ' ', f);
}

static int line(char* s, FILE* f) {
  int n;

  if (fgets(s, N, f) == NULL) return 1;
  n = strlen(s);
  if (n > 0 && s[n - 1] == '\n') s[n - 1] = '\0';
  return 0;
}

static int eq(const char* a, const char* b) {
  return strncmp(a, b, N) == 0;
}
