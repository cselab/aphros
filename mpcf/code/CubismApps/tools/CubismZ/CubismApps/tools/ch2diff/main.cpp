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

#include <float.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

//#include "ParIO.h"
#include <omp.h>
//#define _TRANSPOSE_DATA_

#include "ArgumentParser.h"
#include "Reader_WaveletCompression.h"
#include "Reader_WaveletCompression0.h"

int main(int argc, char **argv)
{
	const double init_t0 = omp_get_wtime();

	/* Initialize MPI */
	MPI_Init(&argc, &argv);

	/* MPI variables */
	MPI_Comm comm  = MPI_COMM_WORLD;
	MPI_Info info  = MPI_INFO_NULL;
	MPI::Intracomm& mycomm = MPI::COMM_WORLD;

	const int mpi_rank = mycomm.Get_rank();
	const int mpi_size = mycomm.Get_size();

	const bool isroot = !mpi_rank;

	ArgumentParser argparser(argc, (const char **)argv);

	if (isroot)
		argparser.loud();
	else
		argparser.mute();

	const string inputfile_name1 = argparser("-simdata1").asString("none");
	const string inputfile_name2 = argparser("-simdata2").asString("none");

//	float threshold = 0.1;	// 0.1%
	const double threshold = (double) argparser("-threshold").asDouble(1);
	printf("threshold = %lf\n", threshold);
//	printf("give threshold =");
//	scanf("%f", &threshold);


	if ((inputfile_name1 == "none")||(inputfile_name2 == "none"))
	{
		printf("usage: %s -simdata1 <filename>  -simdata2 <filename> [-swap] [-wtype <wtype>]\n", argv[0]);
		exit(1);
	}

	const bool swapbytes = argparser.check("-swap");
	const int wtype = argparser("-wtype").asInt(1);

	Reader_WaveletCompressionMPI  myreader1(mycomm, inputfile_name1, swapbytes, wtype);
	Reader_WaveletCompressionMPI0 myreader2(mycomm, inputfile_name2, swapbytes, wtype);

	myreader1.load_file();
	myreader2.load_file();
	const double init_t1 = omp_get_wtime();

	const double t0 = omp_get_wtime(); 

	int dim[3], period[3], reorder;
	int coord[3], id;

	int NBX1 = myreader1.xblocks();
	int NBY1 = myreader1.yblocks();
	int NBZ1 = myreader1.zblocks();
	fprintf(stdout, "[1] I found in total %dx%dx%d blocks.\n", NBX1, NBY1, NBZ1);

	int NBX2 = myreader2.xblocks();
	int NBY2 = myreader2.yblocks();
	int NBZ2 = myreader2.zblocks();
	fprintf(stdout, "[2] I found in total %dx%dx%d blocks.\n", NBX2, NBY2, NBZ2);
	
	if ((NBX1 != NBX2) || (NBZ1 != NBZ2) ||(NBZ1 != NBZ2)) {
		printf("Dimensions differ, exiting..\n");
		exit(1);
	}
	
	static Real targetdata1[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];
	static Real targetdata2[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];

	const int nblocks = NBX1*NBY1*NBZ1;
	const int b_end = ((nblocks + (mpi_size - 1))/ mpi_size) * mpi_size; 

	long n = 0;
	double e_inf = 0;
	double e_1 = 0;
	double e_2 = 0;
	
	double n_inf = 0;
	double n_1 = 0;
	double n_2 = 0;

//	double threshold = 0.1;	// 0.1%

//	printf("give threshold =");
//	scanf("%f", &threshold);

	long over_counter = 0;

	double f1, f2;
	double maxdata = -DBL_MAX;
	double mindata =  DBL_MAX;

	for (int b = mpi_rank; b < b_end; b += mpi_size)	
	{
		int z = b / (NBY1 * NBX1);
		int y = (b / NBX1) % NBY1;
		int x = b % NBX1;

		if (b < nblocks)
		{
#if defined(VERBOSE)
			fprintf(stdout, "loading block( %d, %d, %d )...\n", x, y, z); 
#endif
			double zratio1 = myreader1.load_block2(x, y, z, targetdata1);
			double zratio2 = myreader2.load_block2(x, y, z, targetdata2);

			for (int zb = 0; zb < _BLOCKSIZE_; zb++)
				for (int yb = 0; yb < _BLOCKSIZE_; yb++)
					for (int xb = 0; xb < _BLOCKSIZE_; xb++) {
						f1 = (double) targetdata1[zb][yb][xb];
						f2 = (double) targetdata2[zb][yb][xb];
#if 1
						if (f2 > maxdata) maxdata = f2;
						if (f2 < mindata) mindata = f2;

						double v = fabs(f2);

						if (v > n_inf) {
							n_inf = v;
						}
						n_1 += v;
						n_2 += v*v;

						double err = fabs(f1-f2);

						if (err > e_inf) {
#if VERBOSE
							printf("%15.8f vs %15.8f -> %15.8f (rel_err = %15.8f%%)\n", f1, f2, err, 100.0*(err/v));
#endif
							e_inf = err;
						}
						e_1 += err;
						e_2 += err*err;
						n++;	// number id

						double rel_err;

						if (v != 0.0)
//						if (fabs(v) >= 1e-6)
							rel_err = 100.0 * (err / v);
						else
							rel_err = 0.0;

						// absolute relative % error
						if (rel_err >= threshold)
//						if ((rel_err >= threshold)&&(err >= 0.001)) //&&(v > 0.001))
						{
							over_counter++;
							//printf("%5d: %15.8f %15.8f (%15.8f)\n", over_counter, f1, f2, err); 
							//printf("%5d: %15.8f %15.8f (%15.8f)\n", over_counter, f1, f2, rel_err); 

#if 0
							printf("%5d: %15.8f %15.8f (%15.8f - %10.2f%%) [%4d,%4d,%4d]:(%2d,%2d,%2d)\n",
									over_counter, f1, f2, err, rel_err, x, y, z, xb, yb, zb); 
#endif
						}
#endif
					}
	
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

	// MPI_Reduce for all these
	if (!mpi_rank) {
		MPI_Reduce(MPI_IN_PLACE, &mindata, 	1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
		MPI_Reduce(MPI_IN_PLACE, &maxdata, 	1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(MPI_IN_PLACE, &e_inf, 	1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(MPI_IN_PLACE, &n_inf, 	1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(MPI_IN_PLACE, &e_1,   	1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(MPI_IN_PLACE, &n_1,   	1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(MPI_IN_PLACE, &e_2,   	1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(MPI_IN_PLACE, &n_2,   	1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(MPI_IN_PLACE, &n,     	1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(MPI_IN_PLACE, &over_counter,	1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	} else {
		MPI_Reduce(&mindata, &mindata, 		1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
		MPI_Reduce(&maxdata, &maxdata, 		1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(&e_inf, &e_inf, 		1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(&n_inf, &n_inf, 		1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(&e_1, &e_1,   		1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&n_1, &n_1,   		1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&e_2, &e_2,   		1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&n_2, &n_2,   		1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&n, &n,     			1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&over_counter, &over_counter,1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	}
	

	if (!mpi_rank)
	{
		printf("=========================\n");
		printf("e_inf      = %.3e\n", e_inf);
		printf("n_inf      = %.3e\n", n_inf);
		printf("rel(e_inf) = %.3e\n", e_inf/n_inf);
		printf("\n");

		printf("e_1        = %.3e\n", e_1);
		printf("rel(e_1)   = %.3e\n", e_1/n_1);
		printf("mean(e_1)  = %.3e\n", e_1/n);
		printf("\n");

		printf("e_2        = %.3e\n", sqrt(e_2));
		printf("rel(e_2)   = %.3e\n", sqrt(e_2)/sqrt(n_2));
		printf("mean(e_2)  = %.3e\n", sqrt(e_2)/n);
		printf("\n");

		printf("n          = %ld\n", n);
		printf("counter    = %ld\n", over_counter);
		printf("threshold  = %lf%%\n", threshold);
		printf("percentage = %.3f%%\n", 100.0*(float)over_counter/(float)n);

		double linf = e_inf;
		long nall = n;
		double l1 = e_1 / nall;
		const long double mse = e_2 / nall;
		double l2 = sqrt(e_2) / nall;

		double uncompressed_footprint = sizeof(Real) * nall;
		double compressed_footprint = uncompressed_footprint;

#if 1
		FILE *fpz = fopen(inputfile_name1.c_str(), "rb");
		if (fpz == NULL)
                        printf("fp == NULL for %s\n", argv[3]);
		else
		{
			int fd = fileno(fpz); //if you have a stream (e.g. from fopen), not a file descriptor.
			struct stat buf;
			fstat(fd, &buf);
			compressed_footprint = buf.st_size;
			fclose(fpz);
		}
#endif

		printf("compression-rate: %.2f rel-linf-error: %e rel-mean-error: %e\n",
			uncompressed_footprint / compressed_footprint, linf, l1);

		//https://en.wikipedia.org/wiki/Peak_signal-to-noise_ratio
		printf("mindata = %f\n", mindata);
		printf("maxdata = %f\n", maxdata);

		double psnr;

//		if (normalize)
//			psnr = 10 * (double)log10(1. * 1. / mse);
//		else
            		psnr = 20 * log10((maxdata - mindata) / (2 * sqrt(mse)));


		printf("PSNR-vs-BITRATE: %.04f bps %.04f dB\n", compressed_footprint * 8. / nall, psnr);

	}
	
	/* Close/release resources */
	MPI_Finalize();

	return 0;
}
