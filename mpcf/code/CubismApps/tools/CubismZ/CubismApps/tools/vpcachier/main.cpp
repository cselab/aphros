/*
 *  main.cpp
 *  
 *
 *  Created by Diego Rossinelli on 3/27/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#include <iostream>
#include <string>
#include <mpi.h>
#include <numeric>

#include <ArgumentParser.h>

#include "Reader_WaveletCompression.h"
#include "WaveletTexture3D.h"

int main(int argc, const char **  argv)
{	
	MPI::Init_thread(MPI_THREAD_SERIALIZED);
	
	//create my cartesian communicator
	MPI::Intracomm& mycomm = MPI::COMM_WORLD;
	
	const int myrank = mycomm.Get_rank();
	const int pesize = mycomm.Get_size();
	const bool isroot = !myrank;
	
	ArgumentParser argparser(argc, argv);
	
	if (isroot)
		argparser.loud();
	else
		argparser.mute();
	
	const string pathtovpfile = argparser("-vp").asString("ciccio-bello.vpcache");
	const string pathtosimdata = argparser("-simdata").asString("data.channel0");	
	const double minval = argparser("-min").asDouble(0);
	const double maxval = argparser("-max").asDouble(1);
	const double wavelet_threshold  = argparser("-eps").asDouble(0);
	const bool reading = argparser.check("-read");
	const bool halffloat = argparser.check("-f16");
	const bool swap = argparser.check("-swap");
	const bool swapW = argparser.check("-swapW");
	const int wtype_read = argparser("-wtype_read").asInt(1);
	const int wtype_write = argparser("-wtype_write").asInt(1);

/*
	double minval = 0;
	double maxval = 1;

	char *s;
	s = getenv("MINVAL");
	if (s != NULL) minval = atof(s);
	s = getenv("MAXVAL");
	if (s != NULL) maxval = atof(s);
*/

	//printf("minval = %f, maxval = %f\n", minval, maxval);
	
	//just a mini-test for reading
	if (reading)
	{
		if (isroot)
		{
			WaveletTexture3D_Collection texture_collection(pathtovpfile, wtype_read);
			//texture_collection.set_wtype_read(wtype_read);			

			WaveletTexture3D * texture = new WaveletTexture3D;
			
			for(int iz = 0; iz < texture_collection.get_ztextures(); ++iz)
				for(int iy = 0; iy < texture_collection.get_ytextures(); ++iy)
					for(int ix = 0; ix < texture_collection.get_xtextures(); ++ix)
						texture_collection.read(ix, iy, iz, *texture);
			
			delete texture;
			
			printf("Bella li, tutto in regola. Sa vedum!\n");
		}
	}
	else //ok we are ready to start
	{
		//input is the simulation data
		Reader_WaveletCompressionMPI myreader(mycomm, pathtosimdata, swap, wtype_read);
		myreader.load_file();
		
		//figure out the amount of work per rank
		const int xblocks = myreader.xblocks();
		const int yblocks = myreader.yblocks();
		const int zblocks = myreader.zblocks();
		
		const double gridspacing = 1. / max(max(xblocks, yblocks), zblocks) / _BLOCKSIZE_;
		
		MYASSERT(_BLOCKSIZE_ <= _VOXELS_, "_BLOCKSIZE_ = " << _BLOCKSIZE_ << " is bigger than _VOXELS_ =" << _VOXELS_ << "\n");
		
		const int ghosts1side = 2;
		const int puredata1d = _VOXELS_ - 2 * ghosts1side;
		
		const int xtextures = (xblocks * _BLOCKSIZE_ - 2 * ghosts1side) / puredata1d;
		const int ytextures = (yblocks * _BLOCKSIZE_ - 2 * ghosts1side) / puredata1d;
		const int ztextures = (zblocks * _BLOCKSIZE_ - 2 * ghosts1side) / puredata1d;
		
		const int ntextures = xtextures * ytextures * ztextures;
		
		if (isroot)
		{
			double uncompressed_footprint = 4. * xtextures * ytextures * ztextures * powf(_VOXELS_, 3) / 1024. / 1024.;
			
			printf("I am going to create: %d %d %d textures. Uncompressed memory footprint will be %.2f MB.\n", 
				   xtextures, ytextures, ztextures, uncompressed_footprint);
		}
		
		MYASSERT(xtextures > 0 && ytextures > 0 && ztextures > 0, "NUMBER OF VP BLOCKS IS ZERO!");
		
		const int rankwork = std::max(1, (int)std::ceil(ntextures / (double)pesize));
		const int mystart = std::min(ntextures, rankwork * myrank);
		const int myend = std::min(ntextures, rankwork * (myrank + 1));
		const int mytotalwork = myend - mystart;
		
		//printf("%d) my work is %d - %d\n", myrank, mystart, myend);
		
		//spit a warning in case we starve
		{
			int leastwork = -1;
			mycomm.Reduce(&mytotalwork, &leastwork, 1, MPI_INT, MPI::MIN, 0);
			
			assert(leastwork >= 0 || !isroot);
			
			if (isroot && leastwork == 0)
				printf("WARNING: watchout because some of the ranks have zero work.\n");
		}
		
		//ok, lets collect some stats
		vector<Real> vmaxval, vminval, vavg;
		
		WaveletTexture3D_CollectionMPI texture_collection(mycomm, pathtovpfile, xtextures, ytextures, ztextures, wavelet_threshold, halffloat, swapW, wtype_read, wtype_write);
		
		for(int g = mystart; g < myend; ++g)
		{
			const int gx = g % xtextures;
			const int gy = (g / xtextures) % ytextures;
			const int gz = g / (xtextures * ytextures);
			
			assert(gx >= 0 && gx < xtextures);
			assert(gy >= 0 && gy < ytextures);
			assert(gz >= 0 && gz < ztextures);
			
			WaveletTexture3D * texture = new WaveletTexture3D;
			WaveletsOnInterval::FwtAp * const ptrtexdata = & texture->data()[0][0][0];
			
			const int xdatastart = gx * puredata1d;
			const int ydatastart = gy * puredata1d;
			const int zdatastart = gz * puredata1d;
			
			const int xdataend = (gx + 1) * puredata1d + 2 * ghosts1side;
			const int ydataend = (gy + 1) * puredata1d + 2 * ghosts1side;
			const int zdataend = (gz + 1) * puredata1d + 2 * ghosts1side;
			
			texture->geometry.setup<0>(xdatastart, xdataend, ghosts1side, gridspacing);
			texture->geometry.setup<1>(ydatastart, ydataend, ghosts1side, gridspacing);
			texture->geometry.setup<2>(zdatastart, zdataend, ghosts1side, gridspacing);
			
			assert(xdataend - xdatastart == _VOXELS_);
			assert(ydataend - ydatastart == _VOXELS_);
			assert(zdataend - zdatastart == _VOXELS_);
			
			const int xblockstart = xdatastart / _BLOCKSIZE_;
			const int yblockstart = ydatastart / _BLOCKSIZE_;
			const int zblockstart = zdatastart / _BLOCKSIZE_;
			
			const int xblockend = (xdataend - 1) / _BLOCKSIZE_ + 1;
			const int yblockend = (ydataend - 1) / _BLOCKSIZE_ + 1;
			const int zblockend = (zdataend - 1) / _BLOCKSIZE_ + 1;
			
			for(int bz = zblockstart; bz < zblockend; ++bz)
				for(int by = yblockstart; by < yblockend; ++by)
					for(int bx = xblockstart; bx < xblockend; ++bx)
					{
						Real data[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];
						
						assert(bx >= 0 && bx < xblocks);
						assert(by >= 0 && by < yblocks);
						assert(bz >= 0 && bz < zblocks);
						
						myreader.load_block(bx, by, bz, data);
						
						const int xdststart = std::max(xdatastart, bx * _BLOCKSIZE_) - xdatastart;
						const int ydststart = std::max(ydatastart, by * _BLOCKSIZE_) - ydatastart;
						const int zdststart = std::max(zdatastart, bz * _BLOCKSIZE_) - zdatastart;
						
						const int xdstend = std::min(xdataend, (bx + 1) * _BLOCKSIZE_) - xdatastart;
						const int ydstend = std::min(ydataend, (by + 1) * _BLOCKSIZE_) - ydatastart;
						const int zdstend = std::min(zdataend, (bz + 1) * _BLOCKSIZE_) - zdatastart;
						
						const int xsrcoffset = xdatastart - bx * _BLOCKSIZE_;
						const int ysrcoffset = ydatastart - by * _BLOCKSIZE_;
						const int zsrcoffset = zdatastart - bz * _BLOCKSIZE_;
						
						for(int dz = zdststart; dz < zdstend; ++dz)
							for(int dy = ydststart; dy < ydstend; ++dy)
								for(int dx = xdststart; dx < xdstend; ++dx)
								{
									assert(dx >= 0 && dx < _VOXELS_);
									assert(dy >= 0 && dy < _VOXELS_);
									assert(dz >= 0 && dz < _VOXELS_);
									
									assert(dx + xsrcoffset >= 0 && dx + xsrcoffset < _BLOCKSIZE_);
									assert(dy + ysrcoffset >= 0 && dy + ysrcoffset < _BLOCKSIZE_);
									assert(dz + zsrcoffset >= 0 && dz + zsrcoffset < _BLOCKSIZE_);
									
									ptrtexdata[dx + _VOXELS_ * (dy + _VOXELS_ * dz)] = data[dz + zsrcoffset][dy + ysrcoffset][dx + xsrcoffset];
								}
					}
			
			//normalize the data
			{
				const Real a1 = 1 / (maxval - minval);
				const Real a0 = -minval * a1;
				
				Real mysum = 0, mymax = -HUGE_VAL, mymin = HUGE_VAL;
				
				for(int i = 0; i < _VOXELS_ * _VOXELS_ * _VOXELS_; ++i)
				{
					const float val = ptrtexdata[i]; 
					ptrtexdata[i] = std::min(1.f, std::max(0.f, a0 + a1 * val));
					assert(ptrtexdata[i] >= 0 && ptrtexdata[i] <= 1);
					
					mysum += val;
					mymin = min(mymin, (Real)val);
					mymax = max(mymax, (Real)val);									
				}
				
				vmaxval.push_back(mymax);
				vminval.push_back(mymin);
				vavg.push_back(mysum / (_VOXELS_ * _VOXELS_ * _VOXELS_));
			}
			
			assert(xblockstart >= 0 && xblockstart < xblocks);
			assert(yblockstart >= 0 && yblockstart < yblocks);
			assert(zblockstart >= 0 && zblockstart < zblocks);
			
			assert(xblockend > 0 && xblockend <= xblocks);
			assert(yblockend > 0 && yblockend <= yblocks);
			assert(zblockend > 0 && zblockend <= zblocks);
			
			texture_collection.write(gx, gy, gz, *texture);
			
			delete texture;
		}
		
		//spit some statistics
		{
			//in case we did not do anything, put some fake values in the stat to prevent sigfault
			if (mytotalwork == 0)
			{
				vmaxval.push_back(-HUGE_VAL);
				vminval.push_back(+HUGE_VAL);
				vavg.push_back(0);
			}
			
			float minval = *std::min_element(vminval.begin(), vminval.end());
			float avgval = std::accumulate(vavg.begin(), vavg.end(), 0) / vavg.size();
			float maxval = *std::max_element(vmaxval.begin(), vmaxval.end());
			
			mycomm.Reduce(isroot ? MPI::IN_PLACE : &minval, &minval, 1, MPI_FLOAT, MPI::MIN, 0);
			mycomm.Reduce(isroot ? MPI::IN_PLACE : &avgval, &avgval, 1, MPI_FLOAT, MPI::SUM, 0);
			mycomm.Reduce(isroot ? MPI::IN_PLACE : &maxval, &maxval, 1, MPI_FLOAT, MPI::MAX, 0);
			
			avgval /= mycomm.Get_size();
			
			if (isroot)
				printf("Ok we are done here. The min, avg, max values found are : %.2f %.2f %.2f\n",  minval, avgval, maxval);
		}
		
		if (isroot)
			std::cout << "Also tchuess zaeme gal.\n";
	}
	
	MPI::Finalize();
	
	return 0;
}
