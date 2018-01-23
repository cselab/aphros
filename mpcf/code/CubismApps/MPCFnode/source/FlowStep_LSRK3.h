/*
 *  FlowStep_LSRK3.h
 *  MPCFnode
 *
 *  Created by Diego Rossinelli on 6/15/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <StencilInfo.h>

#include <Diffusion_CPP.h>

#define HPM_Start(x)
#define HPM_Stop(x)

#include "BlockProcessor_OMP.h"
#include "Types.h"

#include "Tests.h"



namespace LSRK3data
{
    extern int verbosity;
    extern int step_id;

    template < typename Kernel, typename Lab >
    struct Diffusion
    {
        StencilInfo stencil;
        Real mu1, mu2, dtinvh;
        int stencil_start[3];
        int stencil_end[3];

        Diffusion(const Real _dtinvh, const Real _mu1=0, const Real _mu2=0): dtinvh(_dtinvh), mu1(_mu1), mu2(_mu2), stencil(-1,-1,-1,2,2,2, true, 6, 0,1,2,3,4,6)
        {
            stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
            stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
        }

        Diffusion(const Diffusion& c): dtinvh(c.dtinvh), mu1(c.mu1), mu2(c.mu2), stencil(-1,-1,-1,2,2,2, true, 6, 0,1,2,3,4,6)
        {
            stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
            stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
        }

        inline void operator()(Lab& lab, const BlockInfo& info, Block_t& o) const
        {
            Kernel kernel(1, mu1, mu2, info.h_gridpoint, 1., dtinvh);

            const Real * const srcfirst = &lab(-1,-1,-1).alpha1rho1;
            const int labSizeRow = lab.template getActualSize<0>();
            const int labSizeSlice = labSizeRow*lab.template getActualSize<1>();
            Real * const destfirst =  &o.tmp[0][0][0][0];
            kernel.compute(srcfirst, Block_t::gptfloats, labSizeRow, labSizeSlice,
                           destfirst, Block_t::gptfloats, Block_t::sizeX, Block_t::sizeX*Block_t::sizeY);
        }
    };

};


class FlowStep_LSRK3
{
private:
    Grid_t& grid;

protected:
    const int verbosity = 1;
    ArgumentParser& parser;

public:
    FlowStep_LSRK3(Grid_t& grid, const Real CFL, ArgumentParser& parser, const int verbosity=1):
    grid(grid), parser(parser), verbosity(verbosity)
    {
        parser.unset_strict_mode();
    }

    //Real operator()(const Real max_dt, const Real current_time=0.0);
};
