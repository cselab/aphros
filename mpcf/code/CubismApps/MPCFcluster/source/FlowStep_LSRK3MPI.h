/*
 *  FlowStep_LSRK3MPI.h
 *  MPCFcluster
 *
 *  Created by Diego Rossinelli on 11/28/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <BlockLabMPI.h>
#include "BlockProcessor_MPI.h"
#include "Types.h"
#include "Tests.h"
#include <StencilInfo.h>
#include <Diffusion_CPP.h>

typedef BlockLabMPI<Lab> LabMPI;
typedef GridMPI< Grid_t > GridMPI_t;

namespace LSRK3data
{
    extern int verbosity;
    extern int step_id;

    template < typename Kernel, typename Lab >
    struct Diffusion
    {
        StencilInfo stencil;
        Real dtinvh;
        int stencil_start[3];
        int stencil_end[3];

        Diffusion(const Real _dtinvh)
          : dtinvh(_dtinvh), stencil(-1,-1,-1,2,2,2, true, 6, 0,1,2,3,4,6)
        {
            stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
            stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
        }

        Diffusion(const Diffusion& c)
          : dtinvh(c.dtinvh), stencil(-1,-1,-1,2,2,2, true, 6, 0,1,2,3,4,6)
        {
            stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
            stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
        }

        inline void operator()(Lab& lab, const BlockInfo& info, Block_t& o) const
        {
            Kernel kernel(dtinvh);

            const Real * const srcfirst = &lab(-1,-1,-1).alpha2;
            const int labSizeRow = lab.template getActualSize<0>();
            const int labSizeSlice = labSizeRow*lab.template getActualSize<1>();
            Real * const destfirst =  &o.tmp[0][0][0][0];
            kernel.compute(srcfirst, Block_t::gptfloats, labSizeRow, labSizeSlice,
                           destfirst, Block_t::gptfloats, Block_t::sizeX, Block_t::sizeX*Block_t::sizeY);
        }
    };

}

template<typename TGrid>
class FlowStep_LSRK3MPI 
{
    TGrid & grid;

    const int verbosity;
    ArgumentParser& parser;

public:
    ~FlowStep_LSRK3MPI() {}

    FlowStep_LSRK3MPI(
        TGrid& grid, ArgumentParser& parser, const int verbosity=1)
      : grid(grid), parser(parser), verbosity(verbosity)
    {
        parser.unset_strict_mode();
    }

    Real operator()(const Real max_dt, const Real current_time=0.0);
};
