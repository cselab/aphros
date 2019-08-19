#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <hdf5.h>
#include <h5.h>

#ifndef H5_HAVE_PARALLEL
#error needs parallel HDF5
#endif

static double buf[65*64*64];
enum {
	X, Y, Z};

int
main(int argc, char **argv)
{
	int status, rank;
	hid_t file;
	MPI_Comm comm;
	char path[] = "p.h5";
	int xlo, ylo, zlo, xs, ys, zs;
	xlo = ylo = zlo = 0;
	xs = ys = zs = 64;
	int start[]	 = {
		zlo, ylo, xlo		};
	int extent[] = {
		zs, ys, xs		};
	int size[] = {
		xs, ys, zs		};
	unsigned int maj, min, rel;

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

	H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
	H5get_libversion(&maj, &min, &rel);
	fprintf(stderr, "hdf5: %d.%d.%d\n", maj, min, rel);

	file = h5_open(comm, path);
	if (file < 0) {
		fprintf(stderr, "%s : h5_fcreate failed\n", path);
		exit(2);
	}
	h5_data(file, size, start, extent, buf);

	double h;
	h = 0.1;
	h5_spacing(file, rank == 0, size, h);

	h5_close(file);
	MPI_Finalize();
	return 0;
}
