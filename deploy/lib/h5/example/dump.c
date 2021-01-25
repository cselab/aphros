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
  USED(argc);
  enum { X, Y, Z };
  enum { I = Z, J = Y, K = X };

  int status, size[3];
  double* buf;
  int i, j, k, m;

  argv++;
  if (*argv == NULL) {
    fprintf(stderr, "needs a path without suffix\n");
    return 2;
  }

  h5_silence();
  status = h5_read_hdf(*argv, size, &buf);
  if (status != 0) {
    fprintf(stderr, "%s:%d: can't read '%s'\n", __FILE__, __LINE__, *argv);
    return 2;
  }
  fprintf(stderr, "size: %d %d %d\n", size[X], size[Y], size[Z]);
  for (m = k = 0; k < size[K]; k++)
    for (j = 0; j < size[J]; j++)
      for (i = 0; i < size[I]; i++)
        printf("%d %d %d %g\n", i, j, k, buf[m++]);
}
