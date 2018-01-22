/*
 *  main.cpp
 *  
 *
 *  Created by Panagiotis Chatzidoukas on 3/28/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#include <iostream>
#include <string>
#include <sstream>
#include <mpi.h>
#include <hdf5.h>
#include <omp.h>
#include "ParIO.h"
//#define _TRANSPOSE_DATA_
#define _COLLECTIVE_IO_

#include "ArgumentParser.h"
#include "Reader_WaveletCompression.h"

	static Real targetdata[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];
	static Real storedata[_BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_];


int main(int argc, const char **argv)
{
	const double init_t0 = omp_get_wtime();

	/* Initialize MPI */
	MPI::Init_thread(MPI_THREAD_SERIALIZED);

	/* MPI variables */
	MPI_Comm comm  = MPI_COMM_WORLD;
	MPI_Info info  = MPI_INFO_NULL;
	MPI::Intracomm& mycomm = MPI::COMM_WORLD;

	const int mpi_rank = mycomm.Get_rank();
	const int mpi_size = mycomm.Get_size();

	const bool isroot = !mpi_rank;

	ArgumentParser argparser(argc, argv);

	if (isroot)
		argparser.loud();
	else
		argparser.mute();

	const string inputfile_name = argparser("-simdata").asString("none");
	const string datfile_name = argparser("-datfile").asString("none");

	if ((inputfile_name == "none")||(datfile_name == "none"))
	{
		printf("usage: %s -simdata <filename>  -datfile <datbasefilename>\n", argv[0]);
		exit(1);
	}

//	size_t dims[4]; 	/* dataset dimensions */
//	size_t	count[4];	/* hyperslab selection parameters */
//	size_t	offset[4];

        const bool swapbytes = argparser.check("-swap");
        const int wtype = argparser("-wtype_read").asInt(1);

#if 1
	Reader_WaveletCompressionMPI myreader(mycomm, inputfile_name, swapbytes, wtype);
#else
	Reader_WaveletCompression myreader(inputfile_name);
#endif
	myreader.load_file();
	const double init_t1 = omp_get_wtime();


	const double t0 = omp_get_wtime(); 

	string datfile_fullname = datfile_name + ".dat";
        MPI_ParIO pio;
        pio.Init((char *)datfile_fullname.c_str(), _BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*sizeof(Real));


//	int dim[3], period[3], reorder;
//	int coord[3], id;

	/* Create a new file collectively and release property list identifier. */
//	FILE *file_id = fopen(datfile_fullname.c_str(), "wb");
//	if (file_id == NULL)
//	{
//		printf("file_id == NULL\n");
//		exit(1);
//	}

	int NBX = myreader.xblocks();
	int NBY = myreader.yblocks();
	int NBZ = myreader.zblocks();
	printf("I found in total %dx%dx%d blocks.\n", NBX, NBY, NBZ);
	
	int NX = NBX*_BLOCKSIZE_;
	int NY = NBY*_BLOCKSIZE_;	
	int NZ = NBZ*_BLOCKSIZE_;

	/* Create the dataspace for the dataset.*/
#if defined(_TRANSPOSE_DATA_)
//	dims[0] = NX;
//	dims[1] = NY;
//	dims[2] = NZ;
//	dims[3] = 1;
#else
//	dims[0] = NZ;
//	dims[1] = NY;
//	dims[2] = NX;
//	dims[3] = 1;
#endif

//	count[0] = _BLOCKSIZE_;
//	count[1] = _BLOCKSIZE_;
//	count[2] = _BLOCKSIZE_;
//	count[3] = 1;

	const int nblocks = NBX*NBY*NBZ;
	const int b_end = ((nblocks + (mpi_size - 1))/ mpi_size) * mpi_size; 

	const Real scalefactor = (Real) argparser("-scalefactor").asDouble(1.0);

	for (int b = mpi_rank; b < b_end; b += mpi_size)	
	{
#if defined(_TRANSPOSE_DATA_)
		int x = b / (NBY * NBZ);
		int y = (b / NBZ) % NBY;
		int z = b % NBZ;
#else
		int z = b / (NBY * NBX);
		int y = (b / NBX) % NBY;
		int x = b % NBX;
#endif

		if (b < nblocks)
		{
//			printf("loading block(%d,%d,%d)\n", x, y, z); 
			myreader.load_block2(x, y, z, targetdata);
	
#if defined(_TRANSPOSE_DATA_)
			for (int xb = 0; xb < _BLOCKSIZE_; xb++)
				for (int yb = 0; yb < _BLOCKSIZE_; yb++)
					for (int zb = 0; zb < _BLOCKSIZE_; zb++)
						storedata[xb*_BLOCKSIZE_*_BLOCKSIZE_ + yb*_BLOCKSIZE_ + zb] = targetdata[zb][yb][xb];
//			offset[0] = x * count[0];
//			offset[1] = y * count[1];
//			offset[2] = z * count[2];
//			offset[3] = 0;
#else
//			offset[0] = z * count[0];
//			offset[1] = y * count[1];
//			offset[2] = x * count[2];
//			offset[3] = 0;
#endif
		

			/* Select hyperslab in the file */
#if defined(_TRANSPOSE_DATA_)
//			fwrite(storedata, sizeof(Real), _BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_, file_id);
			pio.Write((char *)storedata, b);
#else
//			printf("scalefactor = %lf\n" , scalefactor);
			for (int zb = 0; zb < _BLOCKSIZE_; zb++)
				for (int yb = 0; yb < _BLOCKSIZE_; yb++)
					for (int xb = 0; xb < _BLOCKSIZE_; xb++)
						targetdata[zb][yb][xb] /= scalefactor;

//			fwrite(targetdata, sizeof(Real), _BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_, file_id);
			pio.Write((char *)targetdata, b);
#endif
		}
		else {
			/* Dummy write */
		}
	}
	const double t1 = omp_get_wtime(); 

	if (!mpi_rank)
	{
		fprintf(stderr, "Init time = %.3lf seconds\n", init_t1-init_t0);
		fprintf(stderr, "Elapsed time = %.3lf seconds\n", t1-t0);
		fflush(0);
	}
	
	/* Close/release resources */
	pio.Finalize();
	//fclose(file_id);

	MPI::Finalize();

	return 0;
}
