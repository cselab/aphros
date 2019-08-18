#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <hdf5.h>

int
main(int argc, char **argv)
{
    int status, rank;
    hid_t fapl_id, file_id;
    MPI_Comm comm;
    char path[] = "o.h5";

    comm = MPI_COMM_WORLD;
    status = MPI_Init(&argc, &argv);
    if (status != MPI_SUCCESS) {
        fprintf(stderr, "MPI_Init failed\n");
        exit(2);
    }

    status = MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN);
    if (status != MPI_SUCCESS) {
        fprintf(stderr, "MPI_Errhandler_set failed\n");
        exit(2);
    }

    status = MPI_Comm_rank(comm, &rank);
    if (status != MPI_SUCCESS) {
        fprintf(stderr, "MPI_Comm_rank failed\n");
        exit(2);
    }
    if (rank == 0) {
        fapl_id = H5Pcreate(H5P_FILE_ACCESS);
        if (fapl_id < 0) {
            fprintf(stderr, "H5Pcreate failed\n");
            exit(2);
        }
        file_id = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
        if (file_id < 0) {
            fprintf(stderr, "H5Fcreate failed\n");
            exit(2);
        }
        H5Pclose(fapl_id);
    }

    fprintf(stderr, "rank: %d\n", rank);
    MPI_Finalize();
    return 0;
}
