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

  int i, j, k, m, n, siz[3] = {10, 15, 20};
  double *buf, spa, ori[3]= {-0.5, -0.5, -0.5};
  const char *me;

  USED(argc);
  me = argv[0];
  if (*++argv == NULL) {
      fprintf(stderr, "%s: needs a path without suffix\n", me);
      exit(2);
  }

  spa = 0.1;
  n = siz[X] * siz[Y] * siz[Z];
  if ((buf = malloc(n * sizeof(*buf))) == NULL) {
      fprintf(stderr, "%s: fail to malloc (n = %d)\n", me, n);
      exit(2);
  }

  m = 0;
  for (i = 0; i < siz[X]; i++)
    for (j = 0; j < siz[Y]; j++)
      for (k = 0; k < siz[Z]; k++) {
        buf[m++] = i * j * k;
      }

  char* path = *argv;
  char* name = "u";
  if (h5_serial_hdf(path, siz, buf) != 0) {
      fprintf(stderr, "%s: can't wrtie to '%s'\n", me, path);
      exit(2);
  }
  if (h5_xmf(path, name, ori, spa, siz) != 0) {
      fprintf(stderr, "%s can't wrtie to '%s'\n", me, path);
      exit(2);
  }
}
