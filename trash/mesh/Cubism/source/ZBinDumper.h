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

typedef struct _header_serial
{
	long size[8];
} header_serial;

/*
inline size_t ZZcompress(unsigned char *buf, unsigned len, int layout[4], unsigned *max)
inline size_t ZZdecompress(unsigned char * inputbuf, size_t ninputbytes, int layout[4], unsigned char * outputbuf, const size_t maxsize);
*/

using namespace std;

template<typename TGrid, typename Streamer>
void DumpZBin(const TGrid &grid, const int iCounter, const Real t, const string f_name, const string dump_path=".", const bool bDummy=false)
{
	typedef typename TGrid::BlockType B;

    const std::string fullname = f_name + Streamer::EXT;
	char filename[256];
	FILE *file_id;
	int status;

	static const unsigned int NCHANNELS = Streamer::NCHANNELS;
	const unsigned int NX = grid.getBlocksPerDimension(0)*B::sizeX;
	const unsigned int NY = grid.getBlocksPerDimension(1)*B::sizeY;
	const unsigned int NZ = grid.getBlocksPerDimension(2)*B::sizeZ;

    Real memsize = (NX * NY * NZ * sizeof(Real))/(1024.*1024.*1024.);
    cout << "Allocating " << memsize << " GB of BIN data" << endl;
	Real * array_all = new Real[NX * NY * NZ];

	vector<BlockInfo> vInfo_local = grid.getBlocksInfo();

	static const unsigned int sX = 0;
	static const unsigned int sY = 0;
	static const unsigned int sZ = 0;

	static const unsigned int eX = B::sizeX;
	static const unsigned int eY = B::sizeY;
	static const unsigned int eZ = B::sizeZ;

	sprintf(filename, "%s/%s.zbin", dump_path.c_str(), fullname.c_str());

	file_id = fopen(filename, "w");

	header_serial tag;
    fseek(file_id, sizeof(tag), SEEK_SET);
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

        unsigned int max = local_bytes;
        //	int layout[4] = {NCHANNELS, NX, NY, NZ};
        int layout[4] = {NX, NY, NZ, 1};
        long compressed_bytes = ZZcompress((unsigned char *)array_all, local_bytes, layout, &max);	// "in place"

        printf("Writing %ld bytes of Compressed data (cr = %.2f)\n", compressed_bytes, NX*NY*NZ*sizeof(Real)*NCHANNELS*1.0/compressed_bytes);

        tag.size[ichannel] = compressed_bytes;
        size_t wb_data = fwrite(array_all, 1, compressed_bytes, file_id);
    }

    fseek(file_id, 0, SEEK_SET);
    size_t wb_header = fwrite(&tag.size[0], 1, sizeof(tag), file_id);

	status = fclose(file_id);

	delete [] array_all;
}


template<typename TGrid, typename Streamer>
void ReadZBin(TGrid &grid, const string f_name, const string read_path=".")
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

    static const int sX = 0;
    static const int sY = 0;
    static const int sZ = 0;

	const int eX = B::sizeX;
	const int eY = B::sizeY;
	const int eZ = B::sizeZ;

    sprintf(filename, "%s/%s.zbin", read_path.c_str(), f_name.c_str());

    file_id = fopen(filename, "rb");

    long local_count = NX * NY * NZ * 1;
    long local_bytes = local_count * sizeof(Real);

    header_serial tag;
    size_t rb_header = fread(&tag.size[0], 1, sizeof(tag), file_id);

#if DBG
    printf("HEADER(%d):\n", rank);
    for (int i = 0; i < NCHANNELS; i++) {
        printf("channel %d: %ld\n", i, tag.size[i]);
    }
#endif

    for (unsigned int ichannel = 0; ichannel < NCHANNELS; ichannel++)
    {
#if DBG
        printf("compr. size = %ld\n", tag.size[ichannel]); fflush(0);
#endif

        long compressed_bytes = tag.size[ichannel];
#if DBG
        printf("Reading %ld bytes of Compressed data (cr = %.2f)\n", compressed_bytes, local_bytes*1.0/compressed_bytes);
#endif
        unsigned char *tmp = (unsigned char *) malloc(compressed_bytes+4096);

        size_t rb_data = fread(tmp, 1, compressed_bytes, file_id);

        int layout[4] = {NX, NY, NZ, 1};
        size_t decompressed_bytes = ZZdecompress(tmp, compressed_bytes, layout, (unsigned char *)array_all, local_bytes);
        free(tmp);
#if DBG
        printf("size = %ld (%ld)\n", decompressed_bytes, local_bytes); fflush(0);
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

                        Real * const ptr_input = array_all + (gz + NZ * (gy + NY * gx));

                        streamer.operate(*ptr_input, ix, iy, iz, ichannel);	// output -> point
                    }
        }
    } /* ichannel */

    status = fclose(file_id);
    delete [] array_all;
}
