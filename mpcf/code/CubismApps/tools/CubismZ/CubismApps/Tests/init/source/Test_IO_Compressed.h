/*
 *  Test_IO_Compressed.h
 *  CUBISMTests
 *
 *  Created by Panagiotis Hadjidoukas on 4/4/16.
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "Types.h"
#include "SerializerIO_WaveletCompression_MPI_Simple.h"

#include <ArgumentParser.h>
#include <GridMPI.h>

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

		Real minv = 1e8, maxv = -1e8;
		//#pragma omp parallel for
		for(int i=0; i<(int)vInfo.size(); i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;

			for(int iz=0; iz<FluidBlock::sizeZ; iz++)
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
			for(int ix=0; ix<FluidBlock::sizeX; ix++)
			{
				Real p[3];
				info.pos(p, ix, iy, iz);
				Real d = sqrt(pow(p[0]-0.5,2)+pow(p[1]-0.5,2)+pow(p[2]-0.5,2))-0.1;
				b(ix, iy, iz).phi = d; 
				if (d < minv) minv = d;
				if (d > maxv) maxv = d;
			}
		}

		printf("minv = %lf maxv = %lf\n", minv, maxv);
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
			printf("////////////  TEST IO Compressed MPI ///////////////\n");
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

		std::stringstream streamer;
		streamer<<path;
		streamer<<"/";
		streamer<<"output";
		streamer.setf(ios::dec | ios::right);
		streamer.width(5);
		streamer.fill('0');
		streamer<<step_id;

		double threshold = parser("-threshold").asDouble(0.001);
		int wavelet_type = parser("-wtype").asInt(1);

		//edge_width = parser("-width").asInt(0);
		//printf("edge_width = %d\n", edge_width);


		mywaveletdumper.verbose();
		//mywaveletdumper.set_threshold(1e-3);
		mywaveletdumper.set_threshold(threshold);
		mywaveletdumper.set_wtype_write(wavelet_type);
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
