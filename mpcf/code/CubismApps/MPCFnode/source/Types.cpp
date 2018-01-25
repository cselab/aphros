/*
 *  Types.cpp
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 6/21/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#include "Types.h"

Real Simulation_Environment::extent;
Real Simulation_Environment::extents [3];
bool Simulation_Environment::BC_PERIODIC [3];

// Intel compiler does not allow in-class initialization of constant integal
// types.  Still, compilation is possible but linking fails? :(
const int FluidBlock::sizeX = _BLOCKSIZE_;
const int FluidBlock::sizeY = _BLOCKSIZE_;
const int FluidBlock::sizeZ = _BLOCKSIZE_;

FluidElement operator * (const Real a, FluidElement gp)
{
    FluidElement out;

    out.alpha2     = gp.alpha2 * a;

    return out;
}

FluidElement operator + (FluidElement gpa, FluidElement gpb)
{
    FluidElement out;

    out.alpha2     = gpa.alpha2 + gpb.alpha2;

    return out;
}

FluidElement operator - (FluidElement gpa, FluidElement gpb)
{
    FluidElement out;

    out.alpha2     = gpa.alpha2 - gpb.alpha2;

    return out;
}
