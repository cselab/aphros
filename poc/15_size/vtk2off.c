// Created by Sergey Litvinov on 05.10.2019
// Copyright 2019 ETH Zurich

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "arg.h"

enum { N = 1024 };
static char me[] = "ch.vtk2off";

char* argv0;
static int nv, nt;
static float* r;
static int* t;

static void usg(void) {
  fprintf(stderr, "usage: %s < VTK > OFF\n", me);
  exit(0);
}

#define MALLOC(n, p)                                                           \
  do {                                                                         \
    *(p) = malloc(n * sizeof(**(p)));                                          \
    if (*(p) == NULL) {                                                        \
      fprintf(stderr, "%s:%d: alloc failed, n = %d\n", __FILE__, __LINE__, n); \
      exit(2);                                                                 \
    }                                                                          \
  } while (0)

#define FREAD(n, p, f)                                                      \
  do {                                                                      \
    if ((int)fread(p, sizeof(*(p)), (n), (f)) != (n)) {                     \
      fprintf(                                                              \
          stderr, "%s:%d: failt to read, n = %d\n", __FILE__, __LINE__, n); \
      exit(2);                                                              \
    }                                                                       \
  } while (0)

#define FWRITE(n, p, f)                                                      \
  do {                                                                       \
    if ((int)fwrite(p, sizeof(*(p)), (n), (f)) != (n)) {                     \
      fprintf(                                                               \
          stderr, "%s:%d: failt to write, n = %d\n", __FILE__, __LINE__, n); \
      exit(2);                                                               \
    }                                                                        \
  } while (0)

#define SWAP(n, p) swap(n, sizeof(*(p)), p)

static int eq(const char* a, const char* b) {
  return strncmp(a, b, N) == 0;
}

static int line(char* s, FILE* f) {
  int n;

  if (fgets(s, N, f) == NULL) return 1;
  n = strlen(s);
  if (n > 0 && s[n - 1] == '\n') s[n - 1] = '\0';
  return 0;
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

static int read_vtk(void) {
  FILE* f;
  char s[N];
  int i, *a, *b, *t0;

  f = stdin;
  if (line(s, f) != 0) {
    fprintf(stderr, "%s:%d: failt to read\n", __FILE__, __LINE__);
    exit(2);
  }

  if (!eq(s, "# vtk DataFile Version 2.0")) {
    fprintf(stderr, "%s:%d: not a vtk file: '%s'\n", __FILE__, __LINE__, s);
    exit(2);
  }

  line(s, f);
  line(s, f);
  if (!eq(s, "BINARY")) {
    fprintf(stderr, "%s:%d: expect BINARY, got '%s'\n", __FILE__, __LINE__, s);
    exit(2);
  }

  line(s, f);
  if (!eq(s, "DATASET POLYDATA")) {
    fprintf(
        stderr, "%s:%d: expect DATASET POLYDATA, got '%s'\n", __FILE__,
        __LINE__, s);
    exit(2);
  }

  line(s, f);
  if (sscanf(s, "POINTS %d float", &nv) != 1) {
    fprintf(stderr, "%s:%d: expect POINTS, got '%s'\n", __FILE__, __LINE__, s);
    exit(2);
  }
  MALLOC(3 * nv, &r);
  FREAD(3 * nv, r, f);

  while (line(s, f) == 0 && s[0] == '\0')
    ;
  if (sscanf(s, "POLYGONS %d %*d", &nt) != 1) {
    fprintf(
        stderr, "%s:%d: expect POLYGONS, got '%s'\n", __FILE__, __LINE__, s);
    exit(2);
  }
  MALLOC(4 * nt, &t0);
  MALLOC(3 * nt, &t);
  FREAD(4 * nt, t0, f);
  for (i = 0, a = t, b = t0; i < nt; i++) {
    b++;
    *a++ = *b++;
    *a++ = *b++;
    *a++ = *b++;
  }
  swap(3 * nt, sizeof(*t), t);
  free(t0);
  return 0;
}

static int write_off() {
  FILE* f;
  int i, j, k, *t0;

  f = stdout;
  fputs("OFF BINARY\n", f);
  MALLOC(5 * nt, &t0);
  for (i = j = k = 0; i < nt; i++) {
    t0[j++] = 3;
    t0[j++] = t[k++];
    t0[j++] = t[k++];
    t0[j++] = t[k++];
    t0[j++] = 0;
  }
  SWAP(5 * nt, t0);
  i = nv;
  SWAP(1, &i);
  FWRITE(1, &i, f);
  i = nt;
  SWAP(1, &i);
  FWRITE(1, &i, f);
  i = 0;
  SWAP(1, &i);
  FWRITE(1, &i, f);
  FWRITE(3 * nv, r, f);
  FWRITE(5 * nt, t0, f);

  free(t0);
  return 0;
}

int main(int argc, char** argv) {
  ARGBEGIN {
    case 'h':
      usg();
  }
  ARGEND;

  read_vtk();
  write_off();

  free(r);
  free(t);
}
