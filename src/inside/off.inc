// Created by Sergey Litvinov on 31.01.2021
// Copyright 2021 ETH Zurich

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "memory.h"

int off_read(
    FILE* f, int* status, int* pnt, int** ptri, int* pnv, double** pver) {
  enum { Off, Numbers, Ver, Tri, End };
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
  int t;
  int t0;
  int t1;
  int t2;
  int* tri;
  int v;

  state = Off;
  v = t = 0;
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
      case Off:
        if (strncmp(s, "OFF", SIZE)) goto not_off;
        state = Numbers;
        break;
      case Numbers:
        if (sscanf(s, "%d %d %*d", &nv, &nt) != 2) {
          fprintf(
              stderr, "%s:%d: expecting numbers, got '%s'\n", __FILE__,
              __LINE__, s);
          goto err;
        }
        ver = malloc(3 * nv * sizeof(*ver));
        tri = malloc(3 * nt * sizeof(*tri));
        if (tri == NULL) {
          fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
          goto err;
        }
        state = Ver;
        break;
      case Ver:
        if (sscanf(s, "%lf %lf %lf", &x, &y, &z) != 3) {
          fprintf(
              stderr, "%s:%d: expcting vertices got '%s'\n", __FILE__, __LINE__,
              s);
          goto err;
        }
        ver[v++] = x;
        ver[v++] = y;
        ver[v++] = z;
        if (v == 3 * nv) state = Tri;
        break;
      case Tri:
        cnt = sscanf(s, "%d %d %d %d", &npt, &t0, &t1, &t2);
        if (cnt != 4 || npt != 3) {
          fprintf(
              stderr, "%s:%d: expcting triangle got '%s'\n", __FILE__, __LINE__,
              s);
          goto err;
        }
        tri[t++] = t0;
        tri[t++] = t1;
        tri[t++] = t2;
        if (t == 3 * nt) state = End;
        break;
      case End:
        fprintf(stderr, "%s:%d: extra line '%s'\n", __FILE__, __LINE__, s);
        goto err;
        break;
    }
  }
  if (state != End) {
    fprintf(stderr, "%s:%d: off file is not complite\n", __FILE__, __LINE__);
    goto err;
  }
  *pnv = nv;
  *pnt = nt;
  *ptri = tri;
  *pver = ver;
  *status = 0;
  return 0;
err:
  return 1;
not_off:
  *status = 1;
  return 0;
}

int off_write(int nt, const int* tri, int nv, const double* ver, FILE* f) {
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
