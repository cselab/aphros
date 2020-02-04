#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <h5.h>
#include <h5serial.h>
#include "util.h"

#define dimension 3

enum { X, Y, Z };
enum { I = Z, J = Y, K = X };

int main(int argc, char** argv) {
  MPI_Comm comm, cartcomm;
  int status, mpi_rank, mpi_npe, mpi_coords[3], e, nbuf, level;
  int mpi_dims[dimension] = {0, 0, 0};
  int Period[dimension] = {0, 0, 0};
  int ext[3], siz[3], sta[3];
  double L0, X0, Y0, Z0, ori[3], *buf;
  double x, y, z, Delta;
  int i, j, k, m;
  char* path = "ba";

  level = 4; /* input */
  L0 = 1.0;
  X0 = 1;
  Y0 = 2;
  Z0 = 3;

  comm = MPI_COMM_WORLD;
  status = MPI_Init(&argc, &argv);
  if (status != MPI_SUCCESS) {
    WARN(("MPI_Init failed\n"));
    goto err;
  }
  MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN);
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_npe);
  status = MPI_Dims_create(mpi_npe, dimension, mpi_dims);
  if (status != MPI_SUCCESS) {
    WARN(("mpi_dims_create failed"));
    goto err;
  }
  MPI_Cart_create(MPI_COMM_WORLD, dimension, mpi_dims, Period, 0, &cartcomm);
  MPI_Cart_coords(cartcomm, mpi_rank, dimension, mpi_coords);
  WARN(("mpi_dims: %d %d %d", mpi_dims[X], mpi_dims[Y], mpi_dims[Z]));
  WARN(("mpi_dims: %d %d %d", mpi_coords[X], mpi_coords[Y], mpi_coords[Z]));

  Delta = L0 * (1. / (1 << level) / mpi_dims[X]);
  e = 1 << level;
  nbuf = e * e * e;
  buf = malloc(nbuf * sizeof(*buf));
  if (buf == NULL) {
    WARN(("fail to allocate 'nbuf = %d'", nbuf));
    goto err;
  }
  nbuf = e * e * e;
  for (k = m = 0; k < e; k++)
    for (j = 0; j < e; j++)
      for (i = 0; i < e; i++) {
        x = (i + e * mpi_coords[X] + 0.5) * Delta + X0;
        y = (j + e * mpi_coords[Y] + 0.5) * Delta + Y0;
        z = (k + e * mpi_coords[Z] + 0.5) * Delta + Z0;
        buf[m++] = x * y * z;
      }

  ext[I] = ext[J] = ext[K] = e;
  siz[I] = e * mpi_dims[X];
  siz[J] = e * mpi_dims[Y];
  siz[K] = e * mpi_dims[Z];
  sta[I] = e * mpi_coords[X];
  sta[J] = e * mpi_coords[Y];
  sta[K] = e * mpi_coords[Z];
  ori[X] = X0;
  ori[Y] = Y0;
  ori[Z] = Z0;
  status = h5_xmf(path, "u", ori, Delta, siz);
  if (mpi_rank == 0) fprintf(stderr, "%s.h5\n%s.xmf\n", path, path);
  if (status != 0) {
    WARN(("can't write xmf '%s'", path));
    goto err;
  }
  status = h5_hdf(MPI_COMM_WORLD, path, siz, sta, ext, buf);
  if (status != 0) {
    WARN(("can't write hdf '%s'", path));
    goto err;
  }
  free(buf);
  MPI_Finalize();
  return 0;
err:
  MPI_Abort(comm, 2);
}
