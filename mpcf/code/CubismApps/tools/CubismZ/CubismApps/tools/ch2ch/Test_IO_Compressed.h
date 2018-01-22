/*
 *  Test_IO_Compressed.h
 *  CUBISMTests
 *
 *  Created by Babak Hejazialhosseini on 5/16/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "Types_single.h"
#include "SerializerIO_WaveletCompression_MPI_Simple.h"
#include "Reader_WaveletCompression.h"
#include <ArgumentParser.h>
#include <GridMPI.h>
#define VERBOSE 0

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
	string inputfile_name, outputfile_name;
	int myrank;
    
	SerializerIO_WaveletCompression_MPI_SimpleBlocking<G, StreamerGridPointIterative> mywaveletdumper;

	void _ic(FluidGrid& grid)
	{
		vector<BlockInfo> vInfo = grid.getBlocksInfo();
		Real min_u = 1e8, max_u = -1e8;

#if 1
		MPI_Comm comm  = MPI_COMM_WORLD;
		MPI::Intracomm& mycomm = MPI::COMM_WORLD;
		const bool swapbytes = parser.check("-swap");
		//cout << "inputfile_name: "  << inputfile_name << "\n" ;
		const int wtype_read = parser("-wtype_read").asInt(1);


		Reader_WaveletCompressionMPI myreader(mycomm, inputfile_name, swapbytes, wtype_read);

		myreader.load_file();
		int NBX = myreader.xblocks();
		int NBY = myreader.yblocks();
		int NBZ = myreader.zblocks();
		const int nblocks = NBX*NBY*NBZ;
		//fprintf(stdout, "I found in total %dx%dx%d=%d blocks.\n", NBX, NBY, NBZ, nblocks);

		static Real targetdata[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];

		int b = myrank;
		int z = b / (NBY * NBX);
		int y = (b / NBX) % NBY;
		int x = b % NBX;

		#if defined(VERBOSE)
		fprintf(stdout, "loading block( %d, %d, %d )...\n", x, y, z);
		#endif
		double zratio = myreader.load_block2(x, y, z, targetdata);
		#if defined(VERBOSE)
		fprintf(stdout, "compression ratio was %.2lf\n", zratio);
		#endif
#endif

		//#pragma omp parallel for
		for(int i=0; i<(int)vInfo.size(); i++)
		{
			BlockInfo info = vInfo[i];

			printf("loading info.index=[%d,%d,%d]\n", info.index[0], info.index[1], info.index[2]);

			FluidBlock& b = *(FluidBlock*)info.ptrBlock;

			double zratio = myreader.load_block2(info.index[0], info.index[1], info.index[2], targetdata);

    
			for(int iz=0; iz<FluidBlock::sizeZ; iz++)
				for(int iy=0; iy<FluidBlock::sizeY; iy++)
					for(int ix=0; ix<FluidBlock::sizeX; ix++)
					{
						//int g_ix = info.index[0]*FluidBlock::sizeX + ix;
						//int g_iy = info.index[1]*FluidBlock::sizeY + iy;
						//int g_iz = info.index[2]*FluidBlock::sizeZ + iz;
						//printf("B=[%d,%d,%d]:[%d,%d,%d] -> [%d,%d,%d] -> %d\n", info.index[0], info.index[1], info.index[2], ix, iy, iz, g_ix, g_iy, g_iz, g_iz*nx*ny+g_iy*nx+g_ix);

						//Real val = targetdata[ix][iy][iz]; //drand48();
						Real val = targetdata[iz][iy][ix]; //drand48();

						b(ix, iy, iz).u = val;
						if (val < min_u) min_u = val;
						if (val > max_u) max_u = val;
					}
		}
		printf("min_u = %lf max_u = %lf\n", min_u, max_u);
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

		MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

		inputfile_name = parser("-simdata").asString("none");
		outputfile_name = parser("-outdata").asString("none");

		if ((inputfile_name == "none")||(outputfile_name == "none"))
		{
			printf("usage: %s -simdata <filename> -outdata <filename> [-swap] [-wtype_read <type1>] [-wtype_write <type2>]\n", "tool");
			exit(1);
		}


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
		streamer<<outputfile_name;
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
		if (grid!=NULL) {
			delete grid;
			grid = NULL;
		}
	}
};
