/*
 *  Test_IO_Compressed.h
 *  CubismZ
 *
 *  Created by Panagiotis Hadjidoukas 4/4/16.
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <ArgumentParser.h>

#include "Types.h"
#include "SerializerIO_WaveletCompression_MPI_Simple.h"

#include <GridMPI.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

typedef GridMPI< FluidGrid > G;

class Test_IO_Compressed : public Simulation
{
protected:
	
	int BPDX, BPDY, BPDZ;
	int XPESIZE, YPESIZE, ZPESIZE;
	int step_id;
	bool VERBOSITY;
	G * grid;
	ArgumentParser parser;
    
	SerializerIO_WaveletCompression_MPI_SimpleBlocking<G, StreamerGridPointIterative> mywaveletdumper;

	void _ic(FluidGrid& grid)
	{
		vector<BlockInfo> vInfo = grid.getBlocksInfo();

		const string datafname = parser("-fdata").asString("in.dat");
		FILE *fp = fopen(datafname.c_str(), "rb");
		if (fp == NULL)
		{
			printf("fp == NULL for %s\n", datafname.c_str());
			exit(1);
		}

		const Real scalefactor = (Real) parser("-scalefactor").asDouble(1.0);

		long nelems = BPDX*BPDY*BPDZ*_BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_;
		long bytes = nelems*sizeof(Real);
		long nx = BPDX*_BLOCKSIZE_;
		long ny = BPDY*_BLOCKSIZE_;

		Real *dataset = (Real *)malloc(bytes);
		fread(dataset, nelems, sizeof(Real), fp);

		//#pragma omp parallel for
		for(int i=0; i<(int)vInfo.size(); i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;

			//printf(" block (%d): info: (%d,%d,%d)\n", i, info.index[0], info.index[1], info.index[2]);  
			for(int iz=0; iz<FluidBlock::sizeZ; iz++)
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
			for(int ix=0; ix<FluidBlock::sizeX; ix++)
			{
				int g_ix = info.index[0]*FluidBlock::sizeX + ix;
				int g_iy = info.index[1]*FluidBlock::sizeY + iy;
				int g_iz = info.index[2]*FluidBlock::sizeZ + iz;
				//printf("B=[%d,%d,%d]:[%d,%d,%d] -> [%d,%d,%d] -> %d\n", info.index[0], info.index[1], info.index[2], ix, iy, iz, g_ix, g_iy, g_iz, g_iz*nx*ny+g_iy*nx+g_ix);

				Real val;
				val = dataset[g_iz*nx*ny + g_iy*nx + g_ix];
				b(ix, iy, iz).phi = scalefactor * val;
//				if (val < min_u) min_u = val;
//				if (val > max_u) max_u = val;


//				Real p[3];
//				info.pos(p, ix, iy, iz);

//				Real number;
//				number = dataset[];
//				b(ix, iy, iz).phi = scalefactor * number;
			}
		}
		fclose(fp);
	}

	void _setup_mpi_constants(int& xpesize, int& ypesize, int& zpesize)
	{
		xpesize = parser("-xpesize").asInt(1);
		ypesize = parser("-ypesize").asInt(1);
		zpesize = parser("-zpesize").asInt(1);
	}
    
public:
	const bool isroot;
	
	Test_IO_Compressed(const bool isroot, const int argc, const char ** argv) : isroot(isroot) , grid(NULL), parser(argc,argv), step_id(0) { }
    
	void setup()
	{
		parser.unset_strict_mode();

		BPDX = parser("-bpdx").asInt(1);
		BPDY = parser("-bpdy").asInt(1);
		BPDZ = parser("-bpdz").asInt(1);

		_setup_mpi_constants(XPESIZE, YPESIZE, ZPESIZE);

		if (!isroot)
			VERBOSITY = 0;

		if (VERBOSITY)
		{
			printf("////////////////////////////////////////////////////////////\n");
			printf("////////////      TEST IO Compressed MPI     ///////////////\n");
			printf("////////////////////////////////////////////////////////////\n");
		}

		grid = new G(XPESIZE, YPESIZE, ZPESIZE, BPDX, BPDY, BPDZ);

		assert(grid != NULL);

		_ic(*grid);
		vp(*grid, step_id);
	}


	void vp(G& grid, const int step_id)
	{
		if (isroot) cout << "dumping MPI VP ...\n" ;

		const string path = parser("-fpath").asString(".");
		const int wtype_write = parser("-wtype_write").asInt(1);

		std::stringstream streamer;
		streamer<<path;
		streamer<<"/";
		streamer<<"output";
		streamer.setf(ios::dec | ios::right);
		streamer.width(5);
		streamer.fill('0');
		streamer<<step_id;

		double threshold = parser("-threshold").asDouble(1e-5);

		mywaveletdumper.verbose();
                printf("setting threshold to %f\n", threshold);
                mywaveletdumper.set_threshold(threshold);
                mywaveletdumper.set_wtype_write(wtype_write);
		mywaveletdumper.Write<0>(grid, streamer.str());

		if (isroot) cout << "done" << endl;
	}


	void dispose()
	{
		if (grid!=NULL)
		{
			delete grid;
			grid = NULL;
		}
	}
};
