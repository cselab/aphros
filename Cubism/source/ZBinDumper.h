/*
 *  ZBinDumper.h
 *  Cubism
 *
 *  Created by Panos Hadjidoukas on 3/18/14.
 *  Copyright 2014 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <iostream>
#include <vector>
#include <stdio.h>
#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif

#include "BlockInfo.h"
#include "LosslessCompression.h"

/*
inline size_t ZZcompress(unsigned char *buf, unsigned len, int layout[4], unsigned *max)
inline size_t ZZdecompress(unsigned char * inputbuf, size_t ninputbytes, int layout[4], unsigned char * outputbuf, const size_t maxsize);
*/

using namespace std;

template<typename TGrid, typename Streamer>
void DumpZBin(TGrid &grid, const int iCounter, const Real t, const string f_name, const string dump_path=".")
{
	typedef typename TGrid::BlockType B;

	char filename[256];
	FILE *file_id;
	int status;

	static const unsigned int NCHANNELS = Streamer::NCHANNELS;
	const unsigned int NX = grid.getBlocksPerDimension(0)*B::sizeX;
	const unsigned int NY = grid.getBlocksPerDimension(1)*B::sizeY;
	const unsigned int NZ = grid.getBlocksPerDimension(2)*B::sizeZ;

	cout << "Writing ZBIN file\n";
	cout << "Allocating " << (NX * NY * NZ * sizeof(Real) * NCHANNELS)/(1024.*1024.*1024.) << "GB of BIN data";

	Real * array_all = new Real[NX * NY * NZ * NCHANNELS];

	cout << " ...done\n";

	vector<BlockInfo> vInfo_local = grid.getBlocksInfo();

	const unsigned int sX = 0;
	const unsigned int sY = 0;
	const unsigned int sZ = 0;

	const unsigned int eX = B::sizeX;
	const unsigned int eY = B::sizeY;
	const unsigned int eZ = B::sizeZ;

	sprintf(filename, "%s/%s.zbin", dump_path.c_str(), f_name.c_str());

	file_id = fopen(filename, "w");

#pragma omp parallel for
	for(int i=0; i<(int)vInfo_local.size(); i++)
	{
		BlockInfo& info = vInfo_local[i];
		const unsigned int idx[3] = {info.index[0], info.index[1], info.index[2]};
		B & b = *(B*)info.ptrBlock;
		Streamer streamer(b);

                for(unsigned int ix=sX; ix<eX; ix++)
                  for(unsigned int iy=sY; iy<eY; iy++)
                    for(unsigned int iz=sZ; iz<eZ; iz++)
		      {
					Real output[NCHANNELS];
					for(unsigned int i=0; i<NCHANNELS; ++i)
						output[i] = 0;

					streamer.operate(ix, iy, iz, (Real*)output);

					const unsigned int gx = idx[0]*B::sizeX + ix;
					const unsigned int gy = idx[1]*B::sizeY + iy;
					const unsigned int gz = idx[2]*B::sizeZ + iz;

					Real * const ptr = array_all + NCHANNELS*(gz + NZ * (gy + NY * gx));

					for(unsigned int i=0; i<NCHANNELS; ++i)
						ptr[i] = output[i];
				}
	}

	unsigned int max = NX * NY * NZ * sizeof(Real) * NCHANNELS;

#if 0
	size_t wb = fwrite(array_all, (NX * NY * NZ * sizeof(Real) * NCHANNELS), 1, file_id);
#else
	int layout[4] = {NCHANNELS, NX, NY, NZ};
	size_t zb = ZZcompress((unsigned char *)array_all, NX * NY * NZ * sizeof(Real) * NCHANNELS, layout, &max);
	printf("Writing %ld bytes of Compressed data (cr = %.2f)\n", zb, NX*NY*NZ*sizeof(Real)*NCHANNELS*1.0/zb);

	size_t wb1 = fwrite(&zb, 1, sizeof(zb), file_id);
	size_t wb2 = fwrite(array_all, 1, zb, file_id);
#endif

	cout << "writing done\n";

	status = fclose(file_id);

	cout << "closing done\n";
	delete [] array_all;
	cout << "deallocating done\n";
}


template<typename TGrid, typename Streamer>
void ReadZBin(TGrid &grid, const int iCounter, const string f_name, const string read_path=".")
{
	typedef typename TGrid::BlockType B;

	char filename[256];
	int status;
	FILE *file_id;

	const int NX = grid.getBlocksPerDimension(0)*B::sizeX;
	const int NY = grid.getBlocksPerDimension(1)*B::sizeY;
	const int NZ = grid.getBlocksPerDimension(2)*B::sizeZ;
	static const int NCHANNELS = Streamer::NCHANNELS;

	Real * array_all = new Real[NX * NY * NZ * NCHANNELS];

	vector<BlockInfo> vInfo_local = grid.getBlocksInfo();

	const int sX = 0;
	const int sY = 0;
	const int sZ = 0;

	const int eX = B::sizeX;
	const int eY = B::sizeY;
	const int eZ = B::sizeZ;

	sprintf(filename, "%s/%s.zbin", read_path.c_str(), f_name.c_str());

	file_id = fopen(filename, "rb");

#if 0
	size_t rb = fread(array_all, (NX * NY * NZ * sizeof(Real) * NCHANNELS), 1, file_id);
#else
	size_t zb;
	size_t rb1 = fread(&zb, 1, sizeof(zb), file_id);
	printf("Reading %ld bytes of Compressed data (cr = %.2f)\n", zb, NX*NY*NZ*sizeof(Real)*NCHANNELS*1.0/zb);
	unsigned char *tmp = (unsigned char *) malloc(zb);
	size_t rb2 = fread(tmp, 1, zb, file_id);

	int layout[4] = {NCHANNELS, NX, NY, NZ};
 	size_t zwb = ZZdecompress(tmp, zb, layout, (unsigned char *)array_all, NX * NY * NZ * sizeof(Real) * NCHANNELS);
	free(tmp);
#endif

#pragma omp parallel for
	for(int i=0; i<vInfo_local.size(); i++)
	{
		BlockInfo& info = vInfo_local[i];
		const int idx[3] = {info.index[0], info.index[1], info.index[2]};
		B & b = *(B*)info.ptrBlock;
		Streamer streamer(b);

		for(int iz=sZ; iz<eZ; iz++)
			for(int iy=sY; iy<eY; iy++)
				for(int ix=sX; ix<eX; ix++)
				{
					const int gx = idx[0]*B::sizeX + ix;
					const int gy = idx[1]*B::sizeY + iy;
					const int gz = idx[2]*B::sizeZ + iz;

					Real * const ptr_input = array_all + NCHANNELS*(gz + NZ * (gy + NY * gx));

					streamer.operate(ptr_input, ix, iy, iz);
				}
	}

	status = fclose(file_id);

	delete [] array_all;
}
