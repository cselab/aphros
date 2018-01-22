/*
 *  BlockProcessing.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 5/24/09.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <stdlib.h>
#include <vector>
#include <map>
#include <string>
#include <typeinfo>

#include "BlockInfo.h"

#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/pipeline.h"
#include "tbb/concurrent_queue.h"
#include "tbb/cache_aligned_allocator.h"
#include "tbb/task_scheduler_init.h"

using namespace std;
using namespace tbb;

namespace Environment
{
	inline void setup(int threads=-1)
	{
		static tbb::task_scheduler_init * init = NULL;
		
		if (init == NULL)
		{
			const int nthreads = threads==-1? 1 : threads;
			init = new tbb::task_scheduler_init(nthreads);
			//printf("INITIALIZED THREADS=%d (HINT is %d)\n", nthreads, 1);
		}
	}
}

/**
 * Functor to actually perform the operations on the blocks.
 */
template <typename BlockType, typename ProcessingMT>
class BlockProcessingMT_Simple_TBB
{
	ProcessingMT& processing;
	const BlockInfo * ptrInfos;

public:

	BlockProcessingMT_Simple_TBB(const BlockInfo * ptrInfos_, ProcessingMT& processing_):
	processing(processing_), ptrInfos(ptrInfos_){}

	BlockProcessingMT_Simple_TBB(const BlockProcessingMT_Simple_TBB& p):
	processing(p.processing), ptrInfos(p.ptrInfos){}

	//forbidden
	BlockProcessingMT_Simple_TBB& operator=(const BlockProcessingMT_Simple_TBB& p){abort(); return *this;}

	template <typename BlockedRange>
	void operator()(const BlockedRange& r) const
	{
		const int nBlocks = r.end() - r.begin();
		const BlockInfo* v = ptrInfos + r.begin();

		// copy constructor required for ProcessingMT
		ProcessingMT p = processing;
		for(int iB=0; iB<nBlocks; iB++)
		{
			const BlockInfo& info = v[iB];
			BlockType& block = *(BlockType*)(info.ptrBlock);

			// operator()(const BlockInfo&, BlockType&) required for ProcessingMT
			p(info, block);
		}
	}
}; /* BlockProcessingMT_Simple_TBB */

template <typename Grid, typename Lab, typename ProcessingMT, int nSlots>
class BlockProcessingMT_TBB
{
	ProcessingMT& processing;

	const BlockInfo * ptrInfos;

	concurrent_queue<Lab *>& m_availableLabs;
	
	const Real current_time;

public:
	BlockProcessingMT_TBB(concurrent_queue<Lab *>& availableLabs, const BlockInfo * ptrInfos_, Grid& grid, ProcessingMT& processing_, const Real t=0, const bool tensorial=false):
	processing(processing_),
	ptrInfos(ptrInfos_),
	m_availableLabs(availableLabs),
	current_time(t)
	{
		for(int i=0; i<nSlots; i++)
		{
			Lab * lab = NULL;
			int retCode = availableLabs.try_pop(lab);
			assert(retCode == true);
			assert(lab!=NULL);

			// stencil_start and stencil_end required for ProcessingMT
			lab->prepare(grid, processing_.stencil_start, processing_.stencil_end, tensorial);

			availableLabs.push(lab);
		}
	}

	template <typename BlockedRange>
	void operator()(const BlockedRange& r) const
	{
		typedef typename Grid::BlockType BlockType;

		Lab* lab = NULL;
		while(!m_availableLabs.try_pop(lab));
		assert(lab != NULL);

		const int nBlocks = r.end() - r.begin();
		const BlockInfo* v = ptrInfos + r.begin();

		for(int iB=0; iB<nBlocks; iB++)
		{
			const BlockInfo& info = v[iB];
			BlockType& block = *(BlockType*)info.ptrBlock;

			lab->load(info,current_time);

			// operator()(LabType&, const BlockInfo&, BlockType&) required for ProcessingMT
			processing(*lab, info, block);
		}

		m_availableLabs.push(lab);
	}

	BlockProcessingMT_TBB(const BlockProcessingMT_TBB& p):
	processing(p.processing),
	ptrInfos(p.ptrInfos), m_availableLabs(p.m_availableLabs), current_time(p.current_time) {}

private:
	//forbidden
	BlockProcessingMT_TBB& operator=(const BlockProcessingMT_TBB& p){abort(); return *this;}
}; /* BlockProcessingMT_TBB */


/*
 ===============================================================================================================================
 ===============================================================================================================================
 */

/**
 * Process blocks with tbb (using threads).
 * @see BlockProcessing_TBB::process
 */
template <typename BlockType>
class BlockProcessing_TBB
{
	static map<string, vector<void *> >& s_cachedResources;
	static BlockInfo * s_ptrInfos;
	static int s_nBlocks;

protected:
	template <typename Lab>
	static void _getResources(concurrent_queue<Lab*>& resources, const int nSlots)
	{
		//0. checks
		//1. find if the resources are already in the cache
		//2. if not allocate them, store in the cache
		//3. fill the queue

		//0.
		assert(resources.empty());

		//1.
		string key = typeid(Lab).name();
		map<string, vector<void *> >::const_iterator it = s_cachedResources.find(key);
		const bool bCacheMiss = (it == s_cachedResources.end());

		if (bCacheMiss)
		{
			//2.
			typedef cache_aligned_allocator<Lab> lab_allocator;


			vector<void*>& new_resources = s_cachedResources[key];
			new_resources.reserve(nSlots);

			for(int i=0; i<nSlots; i++)
			{
				Lab* t = lab_allocator().allocate(1);
				Lab* lab = new((void*)(t)) Lab();

				//3.
				new_resources.push_back(lab);
				resources.push(lab);
			}
		}
		else //3.
			for(vector<void *>::const_iterator itElem = it->second.begin(); itElem!=it->second.end(); itElem++)
				resources.push((Lab*)*itElem);
	}

	static const BlockInfo * _prepareBlockInfos(vector<BlockInfo>& vInfo)
	{
		typedef cache_aligned_allocator<BlockInfo> info_allocator;

		const int nCurrBlocks = vInfo.size();

		if (s_nBlocks < nCurrBlocks)
		{
			if (s_ptrInfos != NULL)
				info_allocator().deallocate(s_ptrInfos, s_nBlocks);

			s_ptrInfos = info_allocator().allocate(nCurrBlocks);
			s_nBlocks = nCurrBlocks;
		}

		BlockInfo * currInfo = s_ptrInfos;
		for(vector<BlockInfo>::const_iterator it = vInfo.begin(); it != vInfo.end(); it++, currInfo++)
			*currInfo = *it;

		return s_ptrInfos;
	}

public:
	BlockProcessing_TBB()
	{
	}

	virtual ~BlockProcessing_TBB()
	{
	}

	/**
	 * Process blocks in parallel using parallel_for (see tbb-doc) to split up the work.
	 * @param vInfo         Info of all the blocks (to be processed) in the grid
	 *                      (e.g. result of Grid::getBlocksInfo()).
	 * @param c             Collection of all the blocks (to be processed) in the grid
	 *                      (e.g. result of Grid::getBlockCollection()).
	 * @param p             Functor processing the block.
	 *                      See MRAG::Multithreading::DummySimpleBlockFunctor for details.
	 * @param nGranularity  Granularity for the parallel_for (see tbb-doc).
	 *                      Optional: if not set, auto_partitioner (see tbb-doc) will be used.
	 */
	template <typename Processing, typename Grid>
	static void process(vector<BlockInfo>& vInfo,  Processing& p, Grid& grid,
						int nGranularity = -1)
	{
		const bool bAutomatic = nGranularity<0;

		const BlockInfo* infos = _prepareBlockInfos(vInfo);

		BlockProcessingMT_Simple_TBB<BlockType,Processing> body(infos, p);

		if (bAutomatic)
			parallel_for(blocked_range<size_t>(0,vInfo.size()), body,  auto_partitioner());
		else
			parallel_for(blocked_range<size_t>(0,vInfo.size(), nGranularity), body);
	}

	/**
	 * Process blocks in parallel using parallel_for (see tbb-doc) to split up the work.
	 * @param vInfo         Info of all the blocks (to be processed) in the grid
	 *                      (e.g. result of Grid::getBlocksInfo()).
	 * @param c             Collection of all the blocks (to be processed) in the grid
	 *                      (e.g. result of Grid::getBlockCollection()).
	 * @param b             Info on the boundaries of the grid (e.g. result of Grid::getBoundaryInfo())
	 * @param p             Functor processing the block.
	 *                      See MRAG::Multithreading::DummyBlockFunctor for details.
	 * @param nGranularity  Granularity for the parallel_for (see tbb-doc).
	 *                      Optional: if not set, auto_partitioner (see tbb-doc) will be used.
	 */
	template <typename Lab, typename Grid, typename Processing>
	static void process(vector<BlockInfo>& vInfo, Processing& p, Grid& grid,
						const Real t=0, const bool tensorial=false, int nGranularity = -1)
	{
		const int nSlots= (int)(NTHREADS);
		assert(nSlots>0);

		concurrent_queue<Lab *> resources;
		_getResources(resources, nSlots);

		const BlockInfo* infos = _prepareBlockInfos(vInfo);

		BlockProcessingMT_TBB<Grid, Lab, Processing, nSlots> body(resources, infos, grid, p, t, tensorial) ;

		const bool bAutomatic = nGranularity<0;
		if (bAutomatic)
			parallel_for(blocked_range<size_t>(0,vInfo.size()), body,  auto_partitioner());
		else
			parallel_for(blocked_range<size_t>(0,vInfo.size(), nGranularity), body);
	}
}; /* BlockProcessing_TBB */

template <typename BlockType>
map<string, vector<void *> >& BlockProcessing_TBB<BlockType>::s_cachedResources = *(new map<string, vector<void *> >());

template <typename BlockType>
BlockInfo * BlockProcessing_TBB<BlockType>::s_ptrInfos = NULL;

template <typename BlockType>
int BlockProcessing_TBB<BlockType>::s_nBlocks = 0;


/*
 ===============================================================================================================================
 ===============================================================================================================================
 */

