/*
 *  Source_CPP.cpp
 *  MPCFcore
 *
 *  Created by Ursula Rasthofer
 *  Date: January 2016
 *  Copyright 2015 ETH Zurich. All rights reserved.
 *
 */

#include <cassert>
#include <iostream>
#include <cstdlib>
#include <limits>

#include "common.h"
#include "Source_CPP.h"


void Source_CPP::compute (const Real * const srcfirst, const int srcfloats,
        Real * const dstfirst, const int dstfloats) const
{
    assert(srcfloats==dstfloats);

    const int N = _BLOCKSIZE_ * _BLOCKSIZE_ * _BLOCKSIZE_ * srcfloats;

    for(int i=0; i<N; i+=srcfloats)
    {
        const Real rho = srcfirst[i] + srcfirst[i+1];
        const Real invRho = 1.0/rho;
        Real dotu = 0;


        for(int dim=0; dim<3; dim++)
        {
            dstfirst[i+2+dim] = a*dstfirst[i+2+dim] + (rho*g[dim] + f[dim]) * dt;
            dotu += (g[dim] + f[dim]*invRho) * srcfirst[i+2+dim];
        }
        dstfirst[i+5] = a*dstfirst[i+5] + dotu*dt;
    }

}
