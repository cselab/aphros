#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <h5.h>
#include "util.h"

#define dimension 3

enum {X, Y, Z};
enum {I = Z, J = Y, K = X};

int
main(int argc, char **argv)
{
    MPI_Comm comm, cartcomm;
    int status, mpi_rank, mpi_npe, mpi_coords[3], e, nbuf;
    int mpi_dims[dimension] = {0, 0, 0};
    int Period[dimension] = {0, 0, 0};
    int ext[3], siz[3], sta[3];
    double spa, ori[3], *buf;
    comm = MPI_COMM_WORLD;
    status = MPI_Init(&argc, &argv);
    if (status != MPI_SUCCESS) {
	WARN(("MPI_Init failed\n"));
	goto err;
    }
    //MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN);
    MPI_Comm_rank(comm, &mpi_rank);
    MPI_Comm_size(comm, &mpi_npe);

    status = MPI_Dims_create (mpi_npe, dimension, mpi_dims);
    if (status != MPI_SUCCESS) {
	WARN(("mpi_dims_create failed"));
	goto err;
    }
    MPI_Cart_create(MPI_COMM_WORLD, dimension, mpi_dims, Period, 0, &cartcomm);
    MPI_Cart_coords(cartcomm, mpi_rank, dimension, mpi_coords);
    WARN(("mpi_dims: %d %d %d", mpi_dims[X], mpi_dims[Y], mpi_dims[Z]));
    WARN(("mpi_dims: %d %d %d", mpi_coords[X], mpi_coords[Y], mpi_coords[Z]));

    e = 10;
    ext[I] = ext[J] = ext[K] = e;
    siz[I] = e*mpi_dims[X];
    siz[J] = e*mpi_dims[Y];
    siz[K] = e*mpi_dims[Z];

    sta[I] = e*mpi_coords[X];
    sta[J] = e*mpi_coords[Y];
    sta[K] = e*mpi_coords[Z];

    ori[X] = ori[Y] = ori[Z] = 0;
    spa = 0.1;
    nbuf = e * e * e;
    buf = malloc(nbuf*sizeof(*buf));
    if (buf == NULL) {
	WARN(("fail to allocate 'nbuf = %d'", nbuf));
	goto err;
    }

    int i, j, k, m;
    for (k = m = 0; k < e; k++)
	for (j = 0; j < e; j++)
	    for (i = 0; i < e; i++) {
		buf[m++] = 42;
	    }

    char *path = "ba";
    status = h5_xmf(path, "u", ori, spa, siz);
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
