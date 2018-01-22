/*
 *  BlockProcessor_OMP.h
 *  MPCFnode
 *
 *  Created by Fabian Wermelinger 10/13/2016
 *  Copyright 2015 ETH Zurich. All rights reserved.
 *
 */
#ifndef BLOCKPROCESSOR_OMP_H_YNUMIBPC
#define BLOCKPROCESSOR_OMP_H_YNUMIBPC

#include <cstdio>
#include <omp.h>
#include "Types.h"
#include "StencilInfo.h"

template<typename TLab, typename TKernel, typename TGrid>
void processOMP(TKernel kappa, TGrid& grid, const Real t=0.0, const bool record=false)
{
    StencilInfo stencil = kappa.stencil;

    vector<BlockInfo>& myInfo = grid.getBlocksInfo();
    BlockInfo * ary = &myInfo.front();
    const int N = myInfo.size();

    const int NTH = omp_get_max_threads();
    double total_time[NTH];

    static TLab * labs = NULL;

    if (labs == NULL)
    {
#ifndef NDEBUG
        printf("allocating %d labs\n", NTH); fflush(0);
#endif /* NDEBUG */

        labs = new TLab[NTH];
        for(int i = 0; i < NTH; ++i)
            labs[i].prepare(grid, stencil.sx,  stencil.ex,  stencil.sy,  stencil.ey,  stencil.sz,  stencil.ez, stencil.tensorial);
    }

#pragma omp parallel
    {
#ifdef _USE_NUMA_
        const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
        const int mynode = omp_get_thread_num() / cores_per_node;
        numa_run_on_node(mynode);
#endif

        const int tid = omp_get_thread_num();
        total_time[tid] = 0;

        Timer timer;
        TLab& mylab = labs[tid];

#pragma omp for schedule(runtime)
        for(int i=0; i<N; i++)
        {
            //we want to measure the time spent in ghost reconstruction
            timer.start();
            mylab.load(ary[i], t);
            total_time[tid] += timer.stop();

            kappa(mylab, ary[i], *(typename TGrid::BlockType*)ary[i].ptrBlock);
        }

#ifndef NDEBUG
#pragma omp single
        {
            double min_val = total_time[0], max_val = total_time[0], sum = total_time[0];

            for(int i=1; i<NTH; i++)
            {
                min_val = min(min_val, total_time[i]);
                max_val = max(max_val, total_time[i]);
                sum += total_time[i];
            }

            if (record)
                printf("(min,max,avg) of lab.load() is (%5.10e, %5.10e, %5.10e)\n", min_val, max_val, sum/NTH);
        }
#endif /* NDEBUG */
    }
}

#endif /* BLOCKPROCESSOR_OMP_H_YNUMIBPC */
