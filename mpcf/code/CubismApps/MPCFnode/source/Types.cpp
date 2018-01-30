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

