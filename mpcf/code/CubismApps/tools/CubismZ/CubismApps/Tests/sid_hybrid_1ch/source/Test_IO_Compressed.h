/*
 *  Test_IO_Compressed.h
 *  CUBISMTests
 *
 *  Created by Panos Hadjidoukas on 4/4/16.
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
	int mypeindex[3];
	int VarID;
    
	SerializerIO_WaveletCompression_MPI_SimpleBlocking<G, StreamerGridPointIterative> mywaveletdumper;

	void _ic(FluidGrid& grid)
	{
		vector<BlockInfo> vInfo = grid.getBlocksInfo();
		Real min_uvwp = 1e8, max_uvwp = -1e8;

#if 1
		char *filename = (char *)"../data/data.init";

		///////////////////////////////////////

		FILE *fid = fopen(filename,"rb");
		int dims[4];

		fread(dims, sizeof(int), 4, fid);	// read only the first 4 elements from the file

		// Dimensions
		int nx      = dims[0];
		int ny      = dims[1];
		int nz      = dims[2];
		int nvar    = dims[3];

		int nsize   = nx*ny*nz;

		//  The following lines are necessary to move the file pointer to the correct
		// location befory you start reading in the variables

		double dt, time;

		fread(&dt, sizeof(double), 1, fid);     //skipping over 'dt'
		fread(&time, sizeof(double), 1, fid);   //skipping over 'time'

		char varnames[8][nvar];

		for (int var=0; var<nvar; var++) {
			fread(varnames[var], sizeof(char), 8, fid);     //skipping over the stored names
		}

		if (!isroot)
		for (int var=0; var<nvar; var++) {
			printf("var[%d] = ->%s<-\n", var, varnames[var]);
		}

		// There are usually 5 variables (or nvar number of variables) stored in the data
		// file: U, V, W, P, ZMIX. Each fread call makes the file pointer shift to the end
		// of the number of values read. Hence each subsequent fread call will start after
		// where the previous fread call stopped. Add on extra variables if necessary.

		double *UVWP;

		//printf("peindex = [%d,%d,%d]\n", mypeindex[0], mypeindex[1], mypeindex[2]);

#if 0
		UVWP = (double *)malloc(nsize*sizeof(double));
		fread(UVWP, sizeof(double), nsize, fid);
		int base_z = 0;
#else
		int rowsZ = BPDZ*FluidBlock::sizeZ;
		UVWP = (double *)malloc(rowsZ*nx*ny*sizeof(double));	// _BLOCKSIZE_*ZBLOCKS  x-y slices
		int base_z = rowsZ*mypeindex[2];			// offset in the Z direction
		int base_y = 0;						// offset in the Y direction - todo 
		long header = 4*sizeof(int) + 2*sizeof(double) + 8*nvar*sizeof(char);
		long base = VarID*nsize*sizeof(double);
		long offset = nx*ny*base_z*sizeof(double);	// + nx*y*sizeof(double); - todo
		fseek(fid, header + base + offset, SEEK_SET);
		fread(UVWP, sizeof(double), rowsZ*nx*ny, fid);
#endif
		fclose(fid);
#endif

		//#pragma omp parallel for
		for(int i=0; i<(int)vInfo.size(); i++)
		{
			BlockInfo info = vInfo[i];

			printf("info.index=[%d,%d,%d]\n", info.index[0], info.index[1], info.index[2]);

			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
    
			for(int iz=0; iz<FluidBlock::sizeZ; iz++)
				for(int iy=0; iy<FluidBlock::sizeY; iy++)
					for(int ix=0; ix<FluidBlock::sizeX; ix++)
					{
						int g_ix = info.index[0]*FluidBlock::sizeX + ix;
						int g_iy = info.index[1]*FluidBlock::sizeY + iy;
						int g_iz = info.index[2]*FluidBlock::sizeZ + iz;
						//printf("B=[%d,%d,%d]:[%d,%d,%d] -> [%d,%d,%d] -> %d\n", info.index[0], info.index[1], info.index[2], ix, iy, iz, g_ix, g_iy, g_iz, g_iz*nx*ny+g_iy*nx+g_ix);

						Real val;
						int l_iz = g_iz - base_z;			// global index -> local index

						//val = UVWP[g_iz*nx*ny + g_iy*nx + g_ix];
						val = UVWP[l_iz*nx*ny + g_iy*nx + g_ix];
						b(ix, iy, iz).uvwp = val;
						if (val < min_uvwp) min_uvwp = val;
						if (val > max_uvwp) max_uvwp = val;
					}
		}
		printf("min_uvwp = %lf max_uvwp = %lf\n", min_uvwp, max_uvwp);
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

		VarID = parser("-variable").asInt(0);	// 0:U, 1:V, 2:W, 3:P

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

		//int myindex[3];
		grid->peindex(mypeindex);

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

		//edge_width = parser("-width").asInt(0);
		//printf("edge_width = %d\n", edge_width);


		mywaveletdumper.verbose();
		//mywaveletdumper.set_threshold(1e-3);
		mywaveletdumper.set_threshold(threshold);
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
