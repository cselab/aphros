// Created by Sergey Litvinov on 24.12.2019
// Copyright 2019 ETH Zurich

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <h5.h>
#include <h5serial.h>

enum { dimension = 3 };
#define USED(x) \
  if (x)        \
    ;           \
  else {        \
  }

int main(int argc, char** argv) {
  enum { X, Y, Z };
  enum { I = Z, J = Y, K = X };

  MPI_Comm cartcomm, comm = MPI_COMM_WORLD;
  int i, j, k, m, n, mpi_rank, mpi_npe, mpi_coords[dimension];
  int Period[dimension] = {0, 0, 0}, mpi_dims[dimension] = {0, 0, 0},
      reorder = 0;
  double x, y, z, spa, *buf;
  const char *me, *path, *name = "u";
  double ori[dimension] = {-1.1, -1.1, -1.1};

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
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    fprintf(stderr, "%s: MPI_Init failed\n", me);
    exit(2);
  }
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_npe);
  MPI_Dims_create(mpi_npe, dimension, mpi_dims);
  MPI_Cart_create(comm, dimension, mpi_dims, Period, reorder, &cartcomm);
  MPI_Cart_coords(cartcomm, mpi_rank, dimension, mpi_coords);
  int siz[dimension] = {n * mpi_dims[Z], n * mpi_dims[Y], n * mpi_dims[X]};
  int sta[dimension] = {
      n * mpi_coords[Z], n * mpi_coords[Y], n * mpi_coords[X]};
  int ext[dimension] = {n, n, n};
  spa = -2 * ori[X] / siz[X];
  buf = malloc(n * n * n * sizeof(*buf));
  if (buf == NULL) {
    fprintf(stderr, "%s: fail to malloc (n = %d)\n", me, n);
    exit(2);
  }
  for (k = m = 0; k < n; k++)
    for (j = 0; j < n; j++)
      for (i = 0; i < n; i++) {
        x = ori[X] + spa * (i + n * mpi_coords[X] + 0.5);
        y = ori[Y] + spa * (j + n * mpi_coords[Y] + 0.5);
        z = ori[Z] + spa * (k + n * mpi_coords[Z] + 0.5);
        buf[m++] = x * x + y * y + z * z - 1;
      }

  h5_xmf(path, name, ori, spa, siz);
  h5_hdf(MPI_COMM_WORLD, path, siz, sta, ext, buf);
  free(buf);
  MPI_Finalize();
  return 0;
}
