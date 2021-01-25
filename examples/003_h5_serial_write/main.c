// Created by Sergey Litvinov on 24.12.2019
// Copyright 2019 ETH Zurich

#include <stdio.h>
#include <stdlib.h>
#include <h5serial.h>

#define USED(x) \
  if (x)        \
    ;           \
  else {        \
  }

int main(int argc, char** argv) {
  enum { X, Y, Z };
  enum { I = Z, J = Y, K = X };

  int i, j, k, m, n;
  double x, y, z, *buf, spa, ori[3] = {-1.1, -1.1, -1.1};
  const char *me, *path, *name = "u";

  USED(argc);
  me = *argv;
  n = 0;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
      case 'h':
        fprintf(stderr, "%s -n int PREFIX\n", me);
        exit(2);
        break;
      case 'n':
        argv++;
        if (*argv == NULL) {
          fprintf(stderr, "%s: -n needs an argument\n", me);
          exit(2);
        }
        n = atoi(*argv);
        break;
      default:
        fprintf(stderr, "%s: unknow option\n", me);
        exit(2);
    }

  if (n == 0) {
    fprintf(stderr, "%s: -n is not set\n", me);
    exit(2);
  }
  if ((path = *argv) == NULL) {
    fprintf(stderr, "%s: needs a path without suffix\n", me);
    exit(2);
  }

  int siz[] = {n, n, n};
  spa = -2 * ori[X] / siz[X];
  if ((buf = malloc(n * n * n * sizeof(*buf))) == NULL) {
    fprintf(stderr, "%s: fail to malloc (n = %d)\n", me, n);
    exit(2);
  }
  m = 0;
  for (i = 0; i < siz[X]; i++)
    for (j = 0; j < siz[Y]; j++)
      for (k = 0; k < siz[Z]; k++) {
        x = ori[X] + spa * (i + 0.5);
        y = ori[Y] + spa * (j + 0.5);
        z = ori[Z] + spa * (k + 0.5);
        buf[m++] = x * x + y * y + z * z - 1;
      }
  if (h5_serial_hdf(path, siz, buf) != 0) {
    fprintf(stderr, "%s: can't wrtie to '%s'\n", me, path);
    exit(2);
  }
  if (h5_xmf(path, name, ori, spa, siz) != 0) {
    fprintf(stderr, "%s can't wrtie to '%s'\n", me, path);
    exit(2);
  }
  free(buf);
}
