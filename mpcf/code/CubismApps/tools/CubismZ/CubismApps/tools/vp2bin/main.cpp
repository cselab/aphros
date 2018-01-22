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
#include "ParIO.h"
#include <omp.h>
//#define _TRANSPOSE_DATA_

#include "ArgumentParser.h"
#include "Reader_WaveletCompression.h"

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
	const string binfile_name = argparser("-binfile").asString("none");

	if ((inputfile_name == "none")||(binfile_name == "none"))
	{
		printf("usage: %s -simdata <filename>  -binfile <binbasefilename> \n", argv[0]);
		exit(1);
	}

	Reader_WaveletCompressionMPI myreader(mycomm, inputfile_name);

	myreader.load_file();
	const double init_t1 = omp_get_wtime();

	const double t0 = omp_get_wtime(); 

	string binfile_fullname = binfile_name + ".bin";;

	int dim[3], period[3], reorder;
	int coord[3], id;


	/* Create a new file collectively and release property list identifier. */
        MPI_ParIO pio;
	pio.Init((char *)binfile_fullname.c_str(), _BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*sizeof(Real));

	int NBX = myreader.xblocks();
	int NBY = myreader.yblocks();
	int NBZ = myreader.zblocks();
	fprintf(stdout, "I found in total %dx%dx%d blocks.\n", NBX, NBY, NBZ);
	
	static Real targetdata[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];
	static Real storedata[_BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_];

	const int nblocks = NBX*NBY*NBZ;
	const int b_end = ((nblocks + (mpi_size - 1))/ mpi_size) * mpi_size; 

	const Real scalefactor = (Real) argparser("-scalefactor").asDouble(1.0);
	if (isroot) printf("scalefactor = %lf\n" , scalefactor);

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
#if defined(VERBOSE)
			fprintf(stdout, "loading block( %d, %d, %d )...\n", x, y, z); 
#endif
			double zratio = myreader.load_block2(x, y, z, targetdata);
#if defined(VERBOSE)
			fprintf(stdout, "compression ratio was %.2lf\n", zratio); 
#endif
			
			if (scalefactor != 1.0) {
			for (int zb = 0; zb < _BLOCKSIZE_; zb++)
				for (int yb = 0; yb < _BLOCKSIZE_; yb++)
					for (int xb = 0; xb < _BLOCKSIZE_; xb++)
						targetdata[zb][yb][xb] /= scalefactor;
			}


#if defined(_TRANSPOSE_DATA_)
			for (int xb = 0; xb < _BLOCKSIZE_; xb++)
				for (int yb = 0; yb < _BLOCKSIZE_; yb++)
					for (int zb = 0; zb < _BLOCKSIZE_; zb++)
						storedata[xb*_BLOCKSIZE_*_BLOCKSIZE_ + yb*_BLOCKSIZE_ + zb] = targetdata[zb][yb][xb];
#endif


#if defined(_TRANSPOSE_DATA_)
			pio.Write((char *)storedata, b);
#else
			pio.Write((char *)targetdata, b);
#endif
		}
		else {
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	const double t1 = omp_get_wtime(); 

	if (!mpi_rank)
	{
		fprintf(stdout, "Init time = %.3lf seconds\n", init_t1-init_t0);
		fprintf(stdout, "Elapsed time = %.3lf seconds\n", t1-t0);
		fflush(0);
	}
	

	/* Close/release resources */
        pio.Finalize();
	MPI::Finalize();

	return 0;
}
