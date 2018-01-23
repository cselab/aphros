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

template<typename TGrid>
class FlowStep_LSRK3MPI : public FlowStep_LSRK3
{
    TGrid & grid;

    template<typename Kflow, typename Ksource, typename Kdiff, typename Ksurf, typename Kupdate, typename Kupdate_state, typename Ksharp>
    void _process_LSRK3(TGrid& grid, const Real dt, const Real current_time);

    template<typename Kflow, typename Ksource, typename Kdiff, typename Ksurf, typename Kupdate, typename Kupdate_state, typename Ksharp>
    pair<Real, Real> step(TGrid& grid, vector<BlockInfo>& vInfo, const Real a, const Real b, const Real dt, const Real current_time);

public:

    ~FlowStep_LSRK3MPI() {}

    FlowStep_LSRK3MPI(TGrid & grid, const Real CFL, ArgumentParser& parser, const int verbosity):
    FlowStep_LSRK3(grid, CFL, parser, verbosity), grid(grid) {}

    Real operator()(const Real max_dt, const Real current_time=0.0);
};
#endif
