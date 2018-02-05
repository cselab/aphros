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
inline void process(TKernel rhs, TGrid& grid, const Real t=0.0, const bool record=false)
{
    // TKernel=Diffusion
    // TGrid=MPIGrid
    TKernel myrhs = rhs;

    SynchronizerMPI& Synch = grid.sync(myrhs);

    const int nthreads = omp_get_max_threads();

    std::vector<TLab> labs(nthreads);

    for(size_t i = 0; i < labs.size(); ++i)
        labs[i].prepare(grid, Synch);

    MPI_Barrier(grid.getCartComm());

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
            l.load(ar[i], t);
            rhs(l, ar[i], *(typename TGrid::BlockType*)ar[i].ptrBlock);
        }
    }

    MPI_Barrier(grid.getCartComm());
}

#endif /* BLOCKPROCESSOR_MPI_H_XDPWOHW2 */
