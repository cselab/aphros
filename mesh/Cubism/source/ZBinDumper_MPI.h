/*
 *  ZBinDumper_MPI.h
 *  Cubism
 *
 *  Created by Panos Hadjidoukas on 3/20/14.
 *  Copyright 2014 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <cassert>

#define MAX_MPI_PROCS	(16*1024)	// header: 0.25MB* 8

typedef struct _header
{
	long offset[8];	// 1 for single compression, NCHANNELS otherwise
	long size[8];
} header;

#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif

using namespace std;

#include "BlockInfo.h"
#include "LosslessCompression.h"

/*
inline size_t ZZcompress(unsigned char *buf, unsigned len, int layout[4], unsigned *max)
inline size_t ZZdecompress(unsigned char * inputbuf, size_t ninputbytes, int layout[4], unsigned char * outputbuf, const size_t maxsize);
*/

template<typename TGrid, typename Streamer>
void DumpZBin_MPI(TGrid &grid, const int iCounter, const Real t, const string f_name, const string dump_path=".")
{
	typedef typename TGrid::BlockType B;

	int rank, nranks;
	char filename[256];
	MPI_Status status;
	MPI_File file_id;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nranks);

	int coords[3];
	grid.peindex(coords);

	const unsigned int NX = grid.getResidentBlocksPerDimension(0)*B::sizeX;
	const unsigned int NY = grid.getResidentBlocksPerDimension(1)*B::sizeY;
	const unsigned int NZ = grid.getResidentBlocksPerDimension(2)*B::sizeZ;
	static const unsigned int NCHANNELS = Streamer::NCHANNELS;

	if (rank==0)
	{
		cout << "Writing BIN file, NCHANNELS = " << NCHANNELS << "\n";

		//Real memsize = (NX * NY * NZ * NCHANNELS * sizeof(Real))/(1024.*1024.*1024.);
		Real memsize = (NX * NY * NZ * sizeof(Real))/(1024.*1024.*1024.);
		cout << "Allocating " << memsize << "GB of BIN data per rank (" << memsize*nranks << "GB in total)\n";
	}
//	Real * array_all = new Real[NX * NY * NZ * NCHANNELS];
	Real * array_all = new Real[NX * NY * NZ];

	vector<BlockInfo> vInfo_local = grid.getResidentBlocksInfo();

	static const unsigned int sX = 0;
	static const unsigned int sY = 0;
	static const unsigned int sZ = 0;

	static const unsigned int eX = B::sizeX;
	static const unsigned int eY = B::sizeY;
	static const unsigned int eZ = B::sizeZ;

	sprintf(filename, "%s/%s.zbin", dump_path.c_str(), f_name.c_str());

	int rc = MPI_File_open( MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file_id );
	if (rc) {
		printf("Unable to create ZBIN file\n");
		exit(1);
	}


	long previous_offset = 0;
	header tag; 	// moved here

	for (unsigned int ichannel = 0; ichannel < NCHANNELS; ichannel++)
	{


#pragma omp parallel for
	for(unsigned int i=0; i<vInfo_local.size(); i++)
	{
		BlockInfo& info = vInfo_local[i];
		const unsigned int idx[3] = {info.index[0], info.index[1], info.index[2]};
		B & b = *(B*)info.ptrBlock;
		Streamer streamer(b);

		for(unsigned int ix=sX; ix<eX; ix++)
		{
			const unsigned int gx = idx[0]*B::sizeX + ix;
			for(unsigned int iy=sY; iy<eY; iy++)
			{
				const unsigned int gy = idx[1]*B::sizeY + iy;
				for(unsigned int iz=sZ; iz<eZ; iz++)
				{
					const unsigned int gz = idx[2]*B::sizeZ + iz;

					assert((gz + NZ * (gy + NY * gx)) < NX * NY * NZ);

					Real * const ptr = array_all + (gz + NZ * (gy + NY * gx));

					//Real output[NCHANNELS];
					//for(int i=0; i<NCHANNELS; ++i)
					//	output[i] = 0;
					//for(int i=0; i<NCHANNELS; ++i) ptr[i] = output[i];
					//ptr[0] = output[ichannel];

					//streamer.operate(ix, iy, iz, (Real*)output);	// point -> output, todo: add an extra argument for channel

					Real output;
					streamer.operate(ix, iy, iz, &output, ichannel);	// point -> output,
					ptr[0] = output;


				}
			}
		}
	}

//	long local_count = NX * NY * NZ * NCHANNELS;
	long local_count = NX * NY * NZ * 1;
	long local_bytes =  local_count * sizeof(Real);
	long offset; // global offset

	unsigned int max = local_bytes;
//	int layout[4] = {NCHANNELS, NX, NY, NZ};
	int layout[4] = {NX, NY, NZ, 1};
	long compressed_bytes = ZZcompress((unsigned char *)array_all, local_bytes, layout, &max);	// "in place"
#if DBG
	printf("Writing %ld bytes of Compressed data (cr = %.2f)\n", compressed_bytes, local_bytes*1.0/compressed_bytes);
#endif
	MPI_Exscan( &compressed_bytes, &offset, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

	if (rank == 0) offset = 0;

#if DBG
	printf("rank %d, offset = %ld, size = %ld\n", rank, offset, compressed_bytes); fflush(0);
#endif

//	header tag;
	tag.offset[ichannel] = offset + previous_offset;
	tag.size[ichannel] = compressed_bytes;
#if DBG
	printf("rank %d, offset = %ld, size = %ld\n", rank, tag.offset[ichannel], tag.size[ichannel]); fflush(0);
#endif
	previous_offset = (tag.offset[ichannel] + tag.size[ichannel]);
	MPI_Bcast(&previous_offset, 1, MPI_LONG, nranks-1, MPI_COMM_WORLD);

	long base = MAX_MPI_PROCS*sizeof(tag); 	// full Header

	MPI_File_write_at(file_id, base + tag.offset[ichannel], (char *)array_all, tag.size[ichannel], MPI_CHAR, &status);

	}	/* ichannel */

	MPI_File_write_at(file_id, rank*sizeof(tag), &tag, 2*8, MPI_LONG, &status);

	MPI_File_close(&file_id);
	delete [] array_all;
}

template<typename TGrid, typename Streamer>
void ReadZBin_MPI(TGrid &grid, const string f_name, const string dump_path=".")
{
	typedef typename TGrid::BlockType B;

	int rank, nranks;
	char filename[256];
	MPI_Status status;
	MPI_File file_id;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nranks);

	int coords[3];
	grid.peindex(coords);

	const int NX = grid.getResidentBlocksPerDimension(0)*B::sizeX;
	const int NY = grid.getResidentBlocksPerDimension(1)*B::sizeY;
	const int NZ = grid.getResidentBlocksPerDimension(2)*B::sizeZ;
	static const int NCHANNELS = Streamer::NCHANNELS;

	Real * array_all = new Real[NX * NY * NZ * NCHANNELS];

	vector<BlockInfo> vInfo_local = grid.getResidentBlocksInfo();

	static const int sX = 0;
	static const int sY = 0;
	static const int sZ = 0;

	const int eX = B::sizeX;
	const int eY = B::sizeY;
	const int eZ = B::sizeZ;

	sprintf(filename, "%s/%s.zbin", dump_path.c_str(), f_name.c_str());

	int rc = MPI_File_open( MPI_COMM_SELF, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &file_id );
	if (rc) {
		printf("Unable to read ZBIN file\n");
		exit(1);
	}

//	long local_count = NX * NY * NZ * NCHANNELS;
	long local_count = NX * NY * NZ * 1;
	long local_bytes = local_count * sizeof(Real);
	long offset;


	header tag;
	MPI_File_read_at(file_id, rank*sizeof(tag), &tag, 2*8, MPI_LONG, &status);

#if DBG
	printf("HEADER(%d):\n", rank);
	for (int i = 0; i < NCHANNELS; i++) {
		printf("channel %d: %ld %ld\n", i, tag.offset[i], tag.size[i]);
	}
#endif

	for (unsigned int ichannel = 0; ichannel < NCHANNELS; ichannel++)
	{

	//MPI_File_read_at(file_id, (rank*2+0)*sizeof(long), &offset     , 1, MPI_LONG, &status);
	//MPI_File_read_at(file_id, (rank*2+1)*sizeof(long), &compressed_bytes, 1, MPI_LONG, &status);
#if DBG
	printf("rank %d, offset = %ld, compr. size = %ld\n", rank, tag.offset[ichannel], tag.size[ichannel]); fflush(0);
#endif

	long compressed_bytes = tag.size[ichannel];
#if DBG
	printf("Reading %ld bytes of Compressed data (cr = %.2f)\n", compressed_bytes, local_bytes*1.0/compressed_bytes);
#endif
	unsigned char *tmp = (unsigned char *) malloc(compressed_bytes+4096);

	long base = MAX_MPI_PROCS*sizeof(tag);	// Header
	MPI_File_read_at(file_id, base + tag.offset[ichannel], (char *)tmp, compressed_bytes, MPI_CHAR, &status);

//	int layout[4] = {NCHANNELS, NX, NY, NZ};
	int layout[4] = {NX, NY, NZ, 1};
	size_t decompressed_bytes = ZZdecompress(tmp, compressed_bytes, layout, (unsigned char *)array_all, local_bytes);
	free(tmp);
#if DBG
	printf("rank %d, offset = %ld, size = %ld (%ld)\n", rank, offset, decompressed_bytes, local_bytes); fflush(0);
#endif

	#pragma omp parallel for
	for(int i=0; i<vInfo_local.size(); i++)
	{
		BlockInfo& info = vInfo_local[i];
		const int idx[3] = {info.index[0], info.index[1], info.index[2]};
		B & b = *(B*)info.ptrBlock;
		Streamer streamer(b);

                for(int ix=sX; ix<eX; ix++)
		  for(int iy=sY; iy<eY; iy++)
		    for(int iz=sZ; iz<eZ; iz++)
		      {
					const int gx = idx[0]*B::sizeX + ix;
					const int gy = idx[1]*B::sizeY + iy;
					const int gz = idx[2]*B::sizeZ + iz;

					//Real * const ptr_input = array_all + NCHANNELS*(gz + NZ * (gy + NY * gx));
					Real * const ptr_input = array_all + (gz + NZ * (gy + NY * gx));

					streamer.operate(*ptr_input, ix, iy, iz, ichannel);	// output -> point
					//streamer.operate(ptr_input, ix, iy, iz);	// input -> point (all channels)
				}
	}


	} /* ichannel */

	MPI_File_close(&file_id);
	delete [] array_all;
}
