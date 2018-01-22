/*
 *  TestLabs.cpp
 *  MPCFnode
 *
 *  Created by Fabian Wermelinger on 10/20/16.
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef TESTLABS_CPP_RVSNKJ92
#define TESTLABS_CPP_RVSNKJ92

#include "TestLabs.h"

///////////////////////////////////////////////////////////////////////////////
// Cloud cases
///////////////////////////////////////////////////////////////////////////////
double CloudData::rho1 = 0;
double CloudData::rho2 = 0;
double CloudData::p1 = 0;
double CloudData::p2 = 0;

FluidElement CloudDataAcoustic::boundaryElement0;
double CloudDataAcoustic::p_amplitude;
double CloudDataAcoustic::frequency;
double CloudDataAcoustic::t0;
double CloudDataAcoustic::phase0;
double CloudDataAcoustic::sigma;

///////////////////////////////////////////////////////////////////////////////
// SIC cases
///////////////////////////////////////////////////////////////////////////////
FluidElement SICCloudBCData::boundaryElement0;

///////////////////////////////////////////////////////////////////////////////
// Advection cases
///////////////////////////////////////////////////////////////////////////////
std::vector<FluidElement> AdvectionBCData::boundaryElementVector;

#endif /* TESTLABS_CPP_RVSNKJ92 */
