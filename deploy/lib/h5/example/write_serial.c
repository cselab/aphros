#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <h5.h>
#include <h5serial.h>

#define USED(x) \
  if (x)        \
    ;           \
  else {        \
  }

int main(int argc, char** argv) {
  enum { X, Y, Z };
  enum { I = Z, J = Y, K = X };

  int status, siz[3];
  double *buf, spa, ori[3];
  int i, j, k, m, n;

  USED(argc);
  argv++;
  if (*argv == NULL) {
    fprintf(stderr, "needs a path without suffix\n");
    return 2;
  }

  siz[X] = 10;
  siz[Y] = 15;
  siz[Z] = 20;
  ori[X] = -0.5;
  ori[Y] = -0.5;
  ori[Z] = -0.5;
  spa = 0.1;

  n = siz[X] * siz[Y] * siz[Z];
  buf = malloc(n * sizeof(*buf));
  if (!buf) {
    fprintf(stderr, "fail to malloc (n = %d)\n", n);
    return 2;
  }

  m = 0;
  for (i = 0; i < siz[X]; i++)
    for (j = 0; j < siz[Y]; j++)
      for (k = 0; k < siz[Z]; k++) {
        buf[m++] = i * j * k;
      }

  char* path = *argv;
  char* name = "u";
  h5_silence();
  status = h5_serial_hdf(path, siz, buf);
  if (status != 0) {
    fprintf(stderr, "%s:%d: can't wrtie to '%s'\n", __FILE__, __LINE__, *argv);
    return 2;
  }
  status = h5_xmf(path, name, ori, spa, siz);
  if (status != 0) {
    fprintf(stderr, "%s:%d: can't wrtie to '%s'\n", __FILE__, __LINE__, *argv);
    return 2;
  }
}
