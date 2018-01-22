/*
 *  HDF5Dumper_MPI.h
 *  Cubism
 *
 *  Created by Babak Hejazialhosseini on 5/24/09.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <cassert>

#include <hdf5.h>

#ifdef _FLOAT_PRECISION_
#define HDF_REAL H5T_NATIVE_FLOAT
#else
#define HDF_REAL H5T_NATIVE_DOUBLE
#endif

using namespace std;

void ReadHDF5_MPI(void * data, const std::size_t gnx, const std::size_t gny, const std::size_t gnz, const std::size_t lnx, const std::size_t lny, const std::size_t lnz, const std::size_t nchannels, const string f_name, const string dump_path=".")
{
	int myrank;
	char filename[256];
	herr_t status;
	hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	int coords[3] = {myrank,0,0};

	hsize_t count[4] = {lnz, lny, lnx, nchannels};

	hsize_t dims[4] = {gnz, gny, gnx, nchannels};

	hsize_t offset[4] = {
		coords[2]*lnz,
		coords[1]*lny,
		coords[0]*lnx, 0};

	sprintf(filename, "%s/%s.h5", dump_path.c_str(), f_name.c_str());

	H5open();
	fapl_id = H5Pcreate(H5P_FILE_ACCESS);
	status = H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);
	file_id = H5Fopen(filename, H5F_ACC_RDONLY, fapl_id);
	status = H5Pclose(fapl_id);

    dataset_id = H5Dopen2(file_id, "data", H5P_DEFAULT);
	fapl_id = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(fapl_id, H5FD_MPIO_COLLECTIVE);

	fspace_id = H5Dget_space(dataset_id);
	H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
	mspace_id = H5Screate_simple(4, count, NULL);
	status = H5Dread(dataset_id, HDF_REAL, mspace_id, fspace_id, fapl_id, data);

	status = H5Sclose(mspace_id);
	status = H5Sclose(fspace_id);
	status = H5Dclose(dataset_id);
	status = H5Pclose(fapl_id);
	status = H5Fclose(file_id);
	H5close();

       return;
}

