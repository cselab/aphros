/*
 *  BlockProcessor_MPI.h
 *  MPCFcluster
 *
 *  Created by Fabian Wermelinger 10/13/2016
 *  Copyright 2015 ETH Zurich. All rights reserved.
 *
 */
#ifndef BLOCKPROCESSOR_MPI_H_XDPWOHW2
#define BLOCKPROCESSOR_MPI_H_XDPWOHW2

#include <mpi.h>
#include "SynchronizerMPI.h"
#include <omp.h>

template<typename TLab, typename TKernel, typename TGrid>
inline void process(TKernel& rhs, TGrid& grid, const Real t=0.0, const bool record=false)
{
    // TKernel=Diffusion
    // TGrid=MPIGrid

    // Get synchronizer and perform communication
    SynchronizerMPI& Synch = grid.sync(rhs);

    const int nthreads = omp_get_max_threads();

    // One instance of TLab per thread
    std::vector<TLab> labs(nthreads);

    // Prepare labs: initialize boundaries and other private fields,
    // allocate memory for cacheBlock, set reference to grid.
    // CacheBlock keeps data including halos
    for(size_t i = 0; i < labs.size(); ++i)
        labs[i].prepare(grid, Synch);

    MPI_Barrier(grid.getCartComm());

    // Get list of available blocks. 
    // BlockInfo objects point to blocks in grid
    std::vector<BlockInfo> av = Synch.avail();
    const int s = av.size();
    BlockInfo * ar = &av.front();

#pragma omp parallel num_threads(nthreads)
    {
        int tid = omp_get_thread_num();
        TLab& l = labs[tid];

#pragma omp for schedule(dynamic,1)
        for(size_t i=0; i<s; i++)
        {
            // Load data from grid to lab cache (including halos)
            l.load(ar[i], t);

            // Evaluate kernel using 
            // lab cache block.data as src and
            // grid block.tmp as dst
            rhs(l, ar[i], *(typename TGrid::BlockType*)ar[i].ptrBlock);

            // The reason why blocks are copied to lab cache
            // is to have halos available.
            // BlockType normally has only inner cells.
            // Lab gives access to halo nodes with negative indices.
        }
    }

    MPI_Barrier(grid.getCartComm());
}

#endif /* BLOCKPROCESSOR_MPI_H_XDPWOHW2 */
