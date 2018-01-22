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
    
	SerializerIO_WaveletCompression_MPI_SimpleBlocking<G, StreamerGridPointIterative> mywaveletdumper;

	void _ic(FluidGrid& grid)
	{
		vector<BlockInfo> vInfo = grid.getBlocksInfo();
		Real min_u = 1e8, max_u = -1e8;
		Real min_v = 1e8, max_v = -1e8;
		Real min_w = 1e8, max_w = -1e8;
		Real min_p = 1e8, max_p = -1e8;

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

		for (int var=0; var<nvar; var++) {
			printf("var[%d] = ->%s<-\n", var, varnames[var]);
		}

		// There are usually 5 variables (or nvar number of variables) stored in the data
		// file: U, V, W, P, ZMIX. Each fread call makes the file pointer shift to the end
		// of the number of values read. Hence each subsequent fread call will start after
		// where the previous fread call stopped. Add on extra variables if necessary.

		double *dummy;
		double *U, *V, *W, *P;

		//dummy = malloc(nsize*sizeof(double));
		U = (double *)malloc(nsize*sizeof(double));
		V = (double *)malloc(nsize*sizeof(double));
		W = (double *)malloc(nsize*sizeof(double));
		P = (double *)malloc(nsize*sizeof(double));

		// dummy   = fread(fid,nsize,'real*8','ieee-le');    % will produce a column vector
		//U	  = reshape(dummy,nx,ny,nz);% now turning the column vector into a 3D matrix

		//fread(dummy, sizeof(double), nsize, fid);
		fread(U, sizeof(double), nsize, fid);
		fread(V, sizeof(double), nsize, fid);
		fread(W, sizeof(double), nsize, fid);
		fread(P, sizeof(double), nsize, fid);

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

						val = U[g_iz*nx*ny + g_iy*nx + g_ix];
						b(ix, iy, iz).u = val;
						if (val < min_u) min_u = val;
						if (val > max_u) max_u = val;

						val = V[g_iz*nx*ny + g_iy*nx + g_ix];
						b(ix, iy, iz).v = val;
						if (val < min_v) min_v = val;
						if (val > max_v) max_v = val;

						val = W[g_iz*nx*ny + g_iy*nx + g_ix];
						b(ix, iy, iz).w = val;
						if (val < min_w) min_w = val;
						if (val > max_w) max_w = val;

						val = P[g_iz*nx*ny + g_iy*nx + g_ix];
						b(ix, iy, iz).p = val;
						if (val < min_p) min_p = val;
						if (val > max_p) max_p = val;
					}
		}
		printf("min_u = %lf max_u = %lf\n", min_u, max_u);
		printf("min_v = %lf max_v = %lf\n", min_v, max_v);
		printf("min_w = %lf max_w = %lf\n", min_w, max_w);
		printf("min_p = %lf max_p = %lf\n", min_p, max_p);
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
		mywaveletdumper.Write<1>(grid, streamer.str());
		mywaveletdumper.Write<2>(grid, streamer.str());
		mywaveletdumper.Write<3>(grid, streamer.str());

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
