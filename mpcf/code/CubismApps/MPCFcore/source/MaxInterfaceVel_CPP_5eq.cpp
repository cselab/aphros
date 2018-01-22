/*
 *  MaxInterfaceVel_CPP_5eq.cpp
 *  MPCFcore
 *
 *  Created by Ursula Rasthofer
 *  Date: October 2015
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

#include "common.h"
#include "MaxInterfaceVel_CPP_5eq.h"

Real MaxInterfaceVel_CPP_5eq::compute(const Real * const src, const int gptfloats) const
{
    const int N=_BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*gptfloats;
    Real ivel = 0;

    for(int i=0; i<N; i+=gptfloats)
    {
        assert(!isnan(src[i]));
        assert(!isnan(src[i+1]));
        assert(!isnan(src[i+2]));
        assert(!isnan(src[i+3]));
        assert(!isnan(src[i+4]));
        assert(!isnan(src[i+6]));

        const Real a1r1 = max(static_cast<Real>(0.0),src[i]);
        const Real a2r2 = max(static_cast<Real>(0.0),src[i+1]);
        const Real alpha2 = max(static_cast<Real>(0.0), min(src[i+6], static_cast<Real>(1.0)));
        const Real alpha1 = max(static_cast<Real>(0.0), min(static_cast<Real>(1.0)-alpha2, static_cast<Real>(1.0)));

        const Real ru  = src[i+2];
        const Real rv  = src[i+3];
        const Real rw  = src[i+4];

        const Real rInv = static_cast<Real>(1.0)/(a1r1+a2r2);

        const Real normu = (ru*ru + rv*rv + rw*rw)*rInv*rInv;

        ivel = max(ivel, static_cast<Real>(4.0)*alpha2*alpha1*normu);
    }

    return ivel;
}
