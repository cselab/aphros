#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "err.h"
#include "inside.h"
#include "memory.h"

enum { SIZE = 999 };
int stl_read(
  FILE* f, int* status, int* pnt, int** ptri, int* pnv, double** pver) {
  enum { Stl, FacetStart, OuterStart, FacetEnd, OuterEnd, Vertex, End};
  char line[SIZE];
  char* s;
  double* ver;
  double x;
  double y;
  double z;
  int cnt;
  int i;
  int j;
  int npt;
  int nt;
  int nv;
  int state;
  int t0;
  int t1;
  int t2;
  int* tri;
  int v;
  int cap;
  int iver;

  state = Stl;
  v = 0;
  cap = 10;
  if ((ver = malloc(cap * sizeof *ver)) == NULL) {
      fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
      goto err;
  }
  for (;;) {
    if ((s = fgets(line, SIZE, f)) == NULL) break;
    while (isspace(*s))
      s++; /* leading spaces */
    for (i = 0; s[i] != '\0'; i++)
      if (s[i] == '\n' || s[i] == '#') {
        s[i] = '\0';
        break;
      }
    for (i = j = 0; s[i] != '\0'; i++) /* trailing spaces */
      if (!isspace(s[i])) j = i;
    s[j + 1] = '\0';
    if (s[0] == '\0') /* empty line */
      continue;
    switch (state) {
      case Stl:
        if (strncmp(s, "solid", 5)) goto not_stl;
        state = FacetStart;
        break;
      case FacetStart:
          if (!strncmp(s, "endsolid", 8))
              state = End;
          else if (!strncmp(s, "facet normal", 12))
              state = OuterStart;
          else {
              fprintf(stderr, "%s:%d: expecting 'facet normal' or 'endsolid', got '%s'\n", __FILE__, __LINE__);
              goto err;
          }
          break;
      case OuterStart:
          if (strncmp(s, "outer loop", SIZE)) {
              fprintf(stderr, "%s:%d: expecting 'outer loop', got '%s'\n", __FILE__, __LINE__);
              goto err;
          }
          iver = 0;
          state = Vertex;
          break;
      case Vertex:
          if (strncmp(s, "vertex", 6)) {
              fprintf(stderr, "%s:%d: expecting 'vertex', got '%s'\n", __FILE__, __LINE__);
              goto err;
          }
        if (sscanf(s, "vertex %lf %lf %lf", &x, &y, &z) != 3) {
          fprintf(
              stderr, "%s:%d: expcting vertices got '%s'\n", __FILE__, __LINE__,
              s);
          goto err;
        }
        if (v + 2 >= cap) {
            cap *= 2;
            if ((ver = realloc(ver, cap * sizeof *ver)) == NULL) {
                fprintf(stderr, "%s:%d: realloc failed\n", __FILE__, __LINE__);
                goto err;
            }
        }
        ver[v++] = x;
        ver[v++] = y;
        ver[v++] = z;
        iver ++;
        if (iver == 3)
            state = OuterEnd;
        break;
      case OuterEnd:
          if (strncmp(s, "endloop", SIZE)) {
              fprintf(stderr, "%s:%d: expecting 'endloop', got '%s'\n", __FILE__, __LINE__);
              goto err;
          }
          state = FacetEnd;
          break;
      case FacetEnd:
          if (strncmp(s, "endfacet", SIZE)) {
              fprintf(stderr, "%s:%d: expecting 'endfacet', got '%s'\n", __FILE__, __LINE__);
              goto err;
          }
          state = FacetStart;
          break;
      case End:
        fprintf(stderr, "%s:%d: extra line '%s'\n", __FILE__, __LINE__, s);
        goto err;
        break;          
    }
  }
  if (state != End) {
    fprintf(stderr, "%s:%d: stl file is not complite\n", __FILE__, __LINE__);
    goto err;
  }

  nv = v / 3;
  nt = nv / 3;
  if ((tri = malloc(3 * nt * sizeof *tri)) == NULL) {
      fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
      goto err;
  }
  for (i = 0; i < nv; i++)
      tri[i] = i;
  *pnv = nv;
  *pnt = nt;
  *ptri = tri;
  *pver = ver;
  *status = 0;
  return 0;
err:
  return 1;
not_stl:
  *status = 1;
  return 0;
}

int stl_write(int nt, const int* tri, int nv, const double* ver, FILE* f) {
  int i;
  if (fprintf(f, "OFF\n") < 0) {
    fprintf(stderr, "%s:%d: fail to write\n", __FILE__, __LINE__);
    goto err;
  }
  fprintf(f, "%d %d 0\n", nv, nt);
  for (i = 0; i < nv; i++)
    fprintf(
        f, "%.16g %.16g %.16g\n", ver[3 * i], ver[3 * i + 1], ver[3 * i + 2]);
  for (i = 0; i < nt; i++)
    fprintf(f, "3 %d %d %d\n", tri[3 * i], tri[3 * i + 1], tri[3 * i + 2]);
  return 0;
err:
  return 1;
}
