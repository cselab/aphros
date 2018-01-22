// File       : BlockProcessorMPI.h
// Date       : Fri 01 Apr 2016 05:52:01 PM CEST
// Author     : Fabian Wermelinger
// Description: Process all blocks with MPI
// Copyright 2016 ETH Zurich. All Rights Reserved.
#ifndef BLOCKPROCESSORMPI_H_IKFSZWUJ
#define BLOCKPROCESSORMPI_H_IKFSZWUJ

#include <mpi.h>
#include <omp.h>
#include <vector>

using namespace std;

template<typename TLab, typename Operator, typename TGrid>
inline void BlockProcessorMPI(Operator rhs, TGrid& grid, const Real t=0, const bool record=false)
{
    SynchronizerMPI& Synch = grid.sync(rhs);

    vector<BlockInfo> avail0, avail1;

    const int nthreads = omp_get_max_threads();

    TLab * labs = new TLab[nthreads];

    // Setup the static stencil information for this kernel (operator)
    for(int i = 0; i < nthreads; ++i)
        labs[i].prepare(grid, Synch);

    // process inner blocks
    avail0 = Synch.avail_inner();
    const int Ninner = avail0.size();
    BlockInfo * ary0 = &avail0.front();

#pragma omp parallel num_threads(nthreads)
    {
        int tid = omp_get_thread_num();
        TLab& mylab = labs[tid];

#pragma omp for schedule(dynamic,1)
        for(int i=0; i<Ninner; i++)
        {
            mylab.load(ary0[i], t);
            rhs(mylab, ary0[i], *(FluidBlock*)ary0[i].ptrBlock);
        }
    }

    // process remaining blocks
    avail1 = Synch.avail_halo();
    const int Nhalo = avail1.size();
    BlockInfo * ary1 = &avail1.front();

#pragma omp parallel num_threads(nthreads)
    {
        int tid = omp_get_thread_num();
        TLab& mylab = labs[tid];

#pragma omp for schedule(dynamic,1)
        for(int i=0; i<Nhalo; i++)
        {
            mylab.load(ary1[i], t);
            rhs(mylab, ary1[i], *(FluidBlock*)ary1[i].ptrBlock);
        }
    }

    // clean up
    if(labs!=NULL)
    {
        delete[] labs;
        labs=NULL;
    }

    MPI::COMM_WORLD.Barrier();
}

#endif /* BLOCKPROCESSORMPI_H_IKFSZWUJ */
