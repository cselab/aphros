/*
 *  main.cpp
 *  
 *
 *  Created by Diego Rossinelli on 3/27/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */
#include <ArgumentParser.h>
#include "Reader_WaveletCompression.h"

int main(int argc, const char **  argv)
{
	ArgumentParser argparser(argc, argv);
	
	const string pathtosimdata = argparser("-simdata").asString("data.channel0");
	Reader_WaveletCompression myreader(pathtosimdata);
	myreader.load_file();
	printf("I found in total %dx%dx%d blocks.\n", myreader.xblocks(), myreader.yblocks(), myreader.zblocks());
	
	if (true)
	{
		Real targetdata[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];
		
		printf("OK FINAL TEST: THE DATA\n");
		const int xblocks = myreader.xblocks();
		const int yblocks = myreader.yblocks();
		const int zblocks = myreader.zblocks();
		for(int ibz = 0; ibz< zblocks; ++ibz)
			for(int iby = 0; iby< yblocks; ++iby)
				for(int ibx = 0; ibx< xblocks; ++ibx)
				{
					myreader.load_block(ibx, iby, ibz, targetdata);
					for(int iz = 0; iz< _BLOCKSIZE_; ++iz)
						for(int iy = 0; iy< _BLOCKSIZE_; ++iy)
							for(int ix = 0; ix< _BLOCKSIZE_; ++ix)
							{
								assert(!isnan(targetdata[iz][iy][ix]));
								//printf("%d %d %d: %e\n", ix, iy, iz, targetdata[iz][iy][ix]);
							}
				}
	}
	
	return 0;
}
