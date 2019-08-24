#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <h5read.h>
#include <h5.h>

int
main(int argc, char **argv)
{
    int status, size[3];
    double *buf;
    argv++;
    if (*argv == NULL) {
	fprintf(stderr, "needs a path without suffix\n");
	return 2;
    }

    h5_silence();
    status = h5_read_hdf(*argv, size, &buf);
}
