/*
 *  FlowStep_LSRK3MPI.h
 *  MPCFcluster
 *
 *  Created by Diego Rossinelli on 11/28/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#ifndef _FSLSRK3MPI_H
#define _FSLSRK3MPI_H

#include <BlockLabMPI.h>

#include <FlowStep_LSRK3.h>

#include <Histogram.h>
#include <ParIO.h>
#include "BlockProcessor_MPI.h"

typedef BlockLabMPI<Lab> LabMPI;
typedef GridMPI< Grid_t > GridMPI_t;

namespace LSRK3MPIdata
{
    extern double t_fs, t_up;
    extern double t_synch_fs, t_bp_fs;
    extern int counter, GSYNCH, nsynch;

#ifndef _SEQUOIA_
    extern MPI_ParIO_Group hist_group;
    extern MPI_ParIO hist_update, hist_rhs, hist_stepid, hist_nsync;
#endif

    template<typename Kflow, typename Kupdate>
    void notify(double avg_time_rhs, double avg_time_update, const size_t NBLOCKS, const size_t NTIMES);
}

template<typename TGrid>
class FlowStep_LSRK3MPI : public FlowStep_LSRK3
{
    TGrid & grid;

    Real _computeSOS();

    // for interface sharpening
    Real _computeMaxIVel();

    template<typename Kflow, typename Ksource, typename Kdiff, typename Ksurf, typename Kupdate, typename Kupdate_state, typename Ksharp>
    void _process_LSRK3(TGrid& grid, const Real dt, const Real current_time);

    template<typename Kflow, typename Ksource, typename Kdiff, typename Ksurf, typename Kupdate, typename Kupdate_state, typename Ksharp>
    pair<Real, Real> step(TGrid& grid, vector<BlockInfo>& vInfo, const Real a, const Real b, const Real dt, const Real current_time);

public:

    ~FlowStep_LSRK3MPI()
    {
#ifndef _SEQUOIA_
        LSRK3MPIdata::hist_update.Finalize();
        LSRK3MPIdata::hist_rhs.Finalize();
        LSRK3MPIdata::hist_stepid.Finalize();
        LSRK3MPIdata::hist_nsync.Finalize();
#endif
    }

    FlowStep_LSRK3MPI(TGrid & grid, const Real CFL, ArgumentParser& parser, const int verbosity):
    FlowStep_LSRK3(grid, CFL, parser, verbosity), grid(grid)
    {
        if (verbosity) cout << "GSYNCH " << parser("-gsync").asInt(omp_get_max_threads()) << endl;

#ifndef _SEQUOIA_
        static const int pehflag = 0;
        LSRK3MPIdata::hist_group.Init(8, parser("-report").asInt(50), pehflag, grid.getCartComm());

        LSRK3MPIdata::hist_update.Init("hist_UPDATE.bin", &LSRK3MPIdata::hist_group);
        LSRK3MPIdata::hist_rhs.Init("hist_FLOWSTEP.bin", &LSRK3MPIdata::hist_group);
        LSRK3MPIdata::hist_stepid.Init("hist_STEPID.bin", &LSRK3MPIdata::hist_group);
        LSRK3MPIdata::hist_nsync.Init("hist_NSYNCH.bin", &LSRK3MPIdata::hist_group);
#endif
    }

    Real operator()(const Real max_dt, const Real current_time=0.0);
};
#endif
