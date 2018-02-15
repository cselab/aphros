/*
 *  SerializerIO.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 5/28/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once
#include <fstream>

using namespace std;

template<typename GridType, typename Streamer>
class SerializerIO
{
public:
		
	// Typedefs
	typedef typename GridType::BlockType TBlock;
	typedef typename TBlock::ElementType TElement;
	
	// Virtual methods
	virtual void Write(GridType & inputGrid, string fileName, Streamer streamer = Streamer())
	{
		ofstream output(fileName.c_str(),  ios::out);
		
		output << inputGrid;
		
		const vector<BlockInfo> vInfo = inputGrid.getBlocksInfo();
		for(vector<BlockInfo>::const_iterator it = vInfo.begin(); it!= vInfo.end(); ++it)
			((TBlock*)(it->ptrBlock))->template Write<Streamer>(output, streamer);
	}
	
	virtual void Read(GridType & inputGrid, string fileName, Streamer streamer = Streamer())
	{
		ifstream input(fileName.c_str(), ios::in);
		
		input >> inputGrid;
		
		const vector<BlockInfo> vInfo = inputGrid.getBlocksInfo();
		for(vector<BlockInfo>::const_iterator it = vInfo.begin(); it!= vInfo.end(); ++it)
			((TBlock*)(it->ptrBlock))->template Read<Streamer>(input, streamer);
	}
};

