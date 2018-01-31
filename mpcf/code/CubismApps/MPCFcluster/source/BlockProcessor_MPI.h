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
    TKernel myrhs = rhs;

    SynchronizerMPI& Synch = grid.sync(myrhs);

    vector<BlockInfo> avail0, avail1;

    const int nthreads = omp_get_max_threads();

    TLab * labs = new TLab[nthreads];

    for(int i = 0; i < nthreads; ++i)
        labs[i].prepare(grid, Synch);

    static int rounds = -1;
    static int one_less = 1;
    if (rounds == -1)
    {
        char *s = getenv("MYROUNDS");
        if (s != NULL)
            rounds = atoi(s);
        else
            rounds = 0;

        char *s2 = getenv("USEMAXTHREADS");
        if (s2 != NULL)
            one_less = !atoi(s2);
    }

    MPI_Barrier(grid.getCartComm());


    avail0 = Synch.avail_inner();
    const int Ninner = avail0.size();
    BlockInfo * ary0 = &avail0.front();

    int nthreads_first;
    if (one_less)
        nthreads_first = nthreads-1;
    else
        nthreads_first = nthreads;

    if (nthreads_first == 0) nthreads_first = 1;

    int Ninner_first = (nthreads_first)*rounds;
    if (Ninner_first > Ninner) Ninner_first = Ninner;
    int Ninner_rest = Ninner - Ninner_first;

#pragma omp parallel num_threads(nthreads_first)
    {
        int tid = omp_get_thread_num();
        TLab& mylab = labs[tid];

#pragma omp for schedule(dynamic,1)
        for(int i=0; i<Ninner_first; i++)
        {
            mylab.load(ary0[i], t);
            rhs(mylab, ary0[i], *(typename TGrid::BlockType*)ary0[i].ptrBlock);
        }
    }

    avail1 = Synch.avail_halo();
    const int Nhalo = avail1.size();
    BlockInfo * ary1 = &avail1.front();

#pragma omp parallel num_threads(nthreads)
    {
        int tid = omp_get_thread_num();
        TLab& mylab = labs[tid];

#pragma omp for schedule(dynamic,1)
        for(int i=-Ninner_rest; i<Nhalo; i++)
        {
            if (i < 0)
            {
                int ii = i + Ninner;
                mylab.load(ary0[ii], t);
                //rhs(mylab, ary0[ii], *(typename TGrid::BlockType*)ary0[ii].ptrBlock);
            }
            else
            {
                mylab.load(ary1[i], t);
                //rhs(mylab, ary1[i], *(typename TGrid::BlockType*)ary1[i].ptrBlock);
            }
        }
    }

    if(labs!=NULL)
    {
        delete[] labs;
        labs=NULL;
    }

    MPI_Barrier(grid.getCartComm());
}

#endif /* BLOCKPROCESSOR_MPI_H_XDPWOHW2 */
