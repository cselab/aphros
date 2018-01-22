/*
 *  SerializerIO_VP.h
 *  Cubism
 *
 *  Created by Babak Hejazialhosseini  on 6/15/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "Types.h"
#include "Matrix4D.h"

template<typename GridType, typename Streamer>
class SerializerIO_VP
{
	int rank, BPDX, BPDY, BPDZ;
	
public:
	
	SerializerIO_VP(int BPDX, int BPDY, int BPDZ, int rank=0): 
		BPDX(BPDX), BPDY(BPDY), BPDZ(BPDZ), rank(rank) { }
	
	void Write(GridType & inputGrid, string fileName, Streamer streamer = Streamer())
	{
		typedef typename GridType::BlockType TBlock;
		
		assert(rank==0);
		vector<BlockInfo> vInfo = inputGrid.getBlocksInfo();
		
		{
			FILE * file = fopen((fileName+".txt").c_str(), "w");
			assert(file!=NULL);
			fprintf(file, "sizeof(Matrix4D): %d\n",(int)sizeof(Matrix4D<float,true,std::allocator>)); 
			fprintf(file, "Blocks: %d\n", (int)vInfo.size());
			fprintf(file, "Block size: %d %d %d\n", TBlock::sizeX, TBlock::sizeY, TBlock::sizeZ);
			fprintf(file, "Blocks per dim per cpu: %d %d %d\n", BPDX, BPDY, BPDZ);
			
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				//FIXME: 2 times same thing because of VP Reader From Cubism
				fprintf(file, "Block %d: Tree Index: %d %d %d %d %d %d\n", i, info.index[0], info.index[1], info.index[2], info.index[0], info.index[1], info.index[2]); 
			}
			
			fclose(file);
		}
		
		{
			static const int nChannels = streamer.channels;
			Matrix4D<float, true, std::allocator> * matData = new Matrix4D<float,true,std::allocator>(nChannels, TBlock::sizeX, TBlock::sizeY, TBlock::sizeZ);
			
			FILE * file = fopen((fileName+".grid").c_str(), "wb");
			assert(file!=NULL);
			
			
			Real out[nChannels];
			
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				TBlock & block = *(TBlock *)info.ptrBlock;
				
				for(int iz=0; iz<TBlock::sizeZ; iz++)
					for(int iy=0; iy<TBlock::sizeY; iy++)
						for(int ix=0; ix<TBlock::sizeX; ix++)
						{
							streamer.operate(block.data[iz][iy][ix], out);
							
							for(int channel=0; channel<nChannels; channel++)
								matData->Access(channel, ix, iy, iz) = out[channel];
						}
				
				matData->Serialize(file);
			}
			
			fclose(file);
			delete matData;
		}
	}
};

