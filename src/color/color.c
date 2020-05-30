// Created Sergey Litvinov on 26.08.2019
// Copyright 2019 ETH Zurich

#include "color.h"
#include <stdio.h>
#include <stdlib.h>

#define SIZE(a) (int)(sizeof(a) / sizeof(*(a)))
#define MAX_N 42
enum { FREE = -1 };

static int root[MAX_N * MAX_N * MAX_N];

static void make(int v) {
  root[v] = v;
}

static int find(int v) {
  if (v == root[v]) return v;
  return root[v] = find(root[v]);
}

static void union0(int a, int b) {
  a = find(a);
  b = find(b);
  if (a != b) root[b] = a;
}

int COLOR_color(int n, /**/ int* pcnt, int* a) {
  enum { I, J, K };
  int lbl[MAX_N * MAX_N * MAX_N];
  int i, j, k, m, N;
  int u, v, w;
  int idx, jdx;
  int x, cnt;
  int d[][3] = {
      {0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {1, 0, 0},
      {1, 0, 1}, {1, 1, 0}, {1, 1, 1},
  };
  if (n > MAX_N) {
    fprintf(stderr, "%s:%d: n=%d > MAX_N=%d\n", __FILE__, __LINE__, n, MAX_N);
    exit(1);
  }

  N = n * n * n;
  for (i = 0; i < N; i++) {
    make(i);
    lbl[i] = FREE;
  }
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++)
        for (m = 0; m < SIZE(d); m++) {
          if ((u = i + d[m][I]) >= n) continue;
          if ((v = j + d[m][J]) >= n) continue;
          if ((w = k + d[m][K]) >= n) continue;
          idx = i * n * n + j * n + k;
          jdx = u * n * n + v * n + w;
          if (a[idx] == COLOR_EMPTY) continue;
          if (a[jdx] == COLOR_EMPTY) continue;
          union0(idx, jdx);
        }
  cnt = 0;
  for (i = 0; i < N; i++) {
    if (a[i] == COLOR_EMPTY) continue;
    x = find(i);
    if (lbl[x] == FREE) lbl[x] = cnt++;
    a[i] = lbl[x];
  }
  *pcnt = cnt;
  return 0;
}
