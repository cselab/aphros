/*
 *  BlockProcessingMPI.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 11/18/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>

#include "BlockProcessing.h"
#include "SynchronizerMPI.h"

class BlockProcessingMPI
{
	template<typename TGrid, typename Lab, typename Operator>
	struct TBBWorker
	{
		typedef typename TGrid::BlockType BlockType;
		const vector<BlockInfo>& myInfo;
		Operator myrhs;
		TGrid* grid;
		const SynchronizerMPI& synch;
		Real t;
		
		TBBWorker(const vector<BlockInfo>& vInfo, Operator rhs, TGrid& grid, const SynchronizerMPI& synch, const Real t=0): 
		myrhs(rhs), grid(&grid), myInfo(vInfo), synch(synch), t(t) 
		{}
		
		TBBWorker(const TBBWorker& c): myrhs(c.myrhs), grid(c.grid), myInfo(c.myInfo), t(c.t), synch(c.synch){}
		
		void operator()(blocked_range<int> range) const 
		{
			Lab mylab;
			
			mylab.prepare(*grid, synch);
			
			const BlockInfo * ary = &myInfo.front();
			for(int i=range.begin(); i<range.end(); i++)
			{
				mylab.load(ary[i], t);
				myrhs(mylab, ary[i], *(BlockType*)ary[i].ptrBlock);
			}
		}
		
		static void _process(const vector<BlockInfo>& vInfo, Operator rhs, TGrid& grid, const SynchronizerMPI& synch, const Real t=0) {
			tbb::parallel_for(blocked_range<int>(0, vInfo.size()), TBBWorker(vInfo, rhs, grid, synch, t), auto_partitioner() );
		}
		
		static void _process(const vector<BlockInfo>& vInfo, Operator rhs, TGrid& grid, const SynchronizerMPI& synch, const Real t, affinity_partitioner& affinitypart) {
			tbb::parallel_for(blocked_range<int>(0, vInfo.size()), TBBWorker(vInfo, rhs, grid, synch, t), affinitypart);
		}
	};
		
public:
	
	template <typename Processing, typename Grid>
	static void process(vector<BlockInfo>& vInfo,  Processing& p, Grid& grid,
						int nGranularity = -1)
	{
		typedef typename Grid::BlockType BlockType;
		const bool bAutomatic = nGranularity<0;
				
		BlockProcessingMT_Simple_TBB<BlockType,Processing> body(&vInfo.front(), p);
		
		if (bAutomatic)
			parallel_for(blocked_range<size_t>(0,vInfo.size()), body,  auto_partitioner());
		else
			parallel_for(blocked_range<size_t>(0,vInfo.size(), nGranularity), body);
	}
	
	template <typename Lab, typename Grid, typename Processing>
	static void process(const vector<BlockInfo>& vInfo, Processing& p, Grid& grid,
						const Real t=0)
	{
		const SynchronizerMPI& SynchronizerMPI = grid.get_SynchronizerMPI(p);
		TBBWorker<Grid, Lab, Processing>::_process(vInfo, p, grid, SynchronizerMPI, t);
	}
	
	template <typename Lab, typename Grid, typename Processing>
	static void process(const vector<BlockInfo>& vInfo, Processing& p, Grid& grid,
						const Real t, affinity_partitioner& affinitypart)
	{
		const SynchronizerMPI& SynchronizerMPI = grid.get_SynchronizerMPI(p);
		TBBWorker<Grid, Lab, Processing>::_process(vInfo, p, grid, SynchronizerMPI, t, affinitypart);
	}
};
