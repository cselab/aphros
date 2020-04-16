#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <h5.h>
#include <h5serial.h>

static double buf[64 * 64 * 64 * 2];
enum { X, Y, Z };

int main(int argc, char** argv) {
  int status, rank;
  MPI_Comm comm;
  char path[] = "p.000", name[] = "p";
  int xlo, ylo, zlo, xs, ys, zs;
  double origin[3] = {0, 0, 0}, spacing;
  xlo = ylo = zlo = 0;
  xs = ys = zs = 64;
  int start[] = {zlo, ylo, xlo};
  int extent[] = {zs, ys, xs};
  int size[] = {xs, ys, zs};

  comm = MPI_COMM_WORLD;
  status = MPI_Init(&argc, &argv);
  if (status != MPI_SUCCESS) {
    fprintf(stderr, "MPI_Init failed\n");
    exit(2);
  }

  MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN);
  status = MPI_Comm_rank(comm, &rank);
  if (status != MPI_SUCCESS) {
    fprintf(stderr, "MPI_Comm_rank failed\n");
    exit(2);
  }
  fprintf(stderr, "rank: %d\n", rank);

  spacing = 0.1;
  h5_silence();
  h5_hdf(comm, path, size, start, extent, buf);
  h5_xmf(path, name, origin, spacing, size);
  fputs(path, stderr);
  fputc('\n', stderr);
  MPI_Finalize();
  return 0;
}
