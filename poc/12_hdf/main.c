#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <hdf5.h>

#ifndef H5_HAVE_PARALLEL
#error needs parallel HDF5
#endif

static hid_t
h5_fcreate(MPI_Comm comm, const char *path)
{
        hid_t plist, file;
        MPI_Info info;
        int status;

        plist = H5Pcreate(H5P_FILE_ACCESS);
        if (plist < 0) {
                fprintf(stderr, "H5Pcreate failed\n");
                return -1;
        }
        status = MPI_Info_create(&info);
        if (status != MPI_SUCCESS) {
                fprintf(stderr, "MPI_Info_create failed\n");
                return -1;
        }

        H5Pset_fapl_mpio(plist, comm, info);
        file = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
        H5Pclose(plist);
        return file;
}

static herr_t
h5_dwrite(hid_t file, hid_t memspace, hid_t filespace, double *buf)
{
        hid_t plist, dset;
        herr_t err;

        dset = H5Dcreate(file, "data", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        plist = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
        err = H5Dwrite(dset, H5T_NATIVE_DOUBLE, memspace, filespace, plist, buf);
        H5Pclose(plist);
        H5Dclose(dset);
        return err;
}

static int
f(hid_t file, int *isize, int *istart, int *iextent, double *buf)
{
        hid_t filespace, memspace;
        herr_t err;
        int dim;
        hsize_t size[] = {
                isize[0], isize[1], isize[2], 1	};
        hsize_t start[] = {
                istart[0], istart[1], istart[2], 0	};
        hsize_t extent[] = {
                iextent[0], iextent[1], iextent[2], 0	};

        dim = 4;
        filespace = H5Screate_simple(dim, size, NULL);
        memspace = H5Screate_simple(dim, extent, NULL);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, extent, NULL);
        err = h5_dwrite(file, memspace, filespace, buf);
        H5Sclose(memspace);
        H5Sclose(filespace);
        return err;
}

static int
g(hid_t file, const char *name, int isize, double *buf)
{
  hid_t space, dset;
        herr_t err;
        int dim;
        hsize_t size[] = {isize};
        dim = 1;
        space = H5Screate_simple(dim, size, NULL);
        dset = H5Dcreate(file, name, H5T_NATIVE_DOUBLE,  space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        err = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,  H5P_DEFAULT, buf);
        H5Sclose(space);
        H5Dclose(dset);
        return err;
}


int
main(int argc, char **argv)
{
        int status, rank;
        hid_t file;
        MPI_Comm comm;
        char path[] = "o.h5";
        int xlo, ylo, zlo, xs, ys, zs;
        xlo = ylo = zlo = 0;
        xs = ys = zs = 10;
        int start[]      = {
                zlo, ylo, xlo};
        int extent[] = {
                zs, ys, xs};
        int size[] = {
                10*xs, 10*ys, 10*zs};
        double buf[10*10*10];
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

        file = h5_fcreate(comm, path);
        if (file < 0) {
                fprintf(stderr, "%s : h5_fcreate failed\n", path);
                exit(2);
        }
        f(file, size, start, extent, buf);
        if (rank == 0) {
          g(file, "vx", 10, buf);
        }
        H5Fclose(file);

        MPI_Finalize();
        return 0;
}
