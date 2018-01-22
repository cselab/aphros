/*
 *  Update_State_CPP_5eq.cpp
 *  MPCFcore
 *
 *  Created by Fabian Wermelinger 03/19/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#include <cassert>
#include <iostream>
#include <cstdlib>
#include <limits>

#include "Update_State_CPP_5eq.h"

void Update_State_CPP_5eq::compute(Real * const dst, const int gptfloats) const
{
    const int N = _BLOCKSIZE_ * _BLOCKSIZE_ * _BLOCKSIZE_ * gptfloats;

    for(int i=0; i<N; i+=gptfloats)
    {
        const Real r1 = dst[i];
        const Real r2 = dst[i+1];
        const Real u  = dst[i+2];
        const Real v  = dst[i+3];
        const Real w  = dst[i+4];
        const Real e  = dst[i+5];
        const Real A2 = dst[i+6];

        assert(!isinf(dst[i]));   assert(!isnan(dst[i]));
        assert(!isinf(dst[i+1])); assert(!isnan(dst[i+1]));
        assert(!isinf(dst[i+2])); assert(!isnan(dst[i+2]));
        assert(!isinf(dst[i+3])); assert(!isnan(dst[i+3]));
        assert(!isinf(dst[i+4])); assert(!isnan(dst[i+4]));
        assert(!isinf(dst[i+5])); assert(!isnan(dst[i+5]));
        assert(!isinf(dst[i+6])); assert(!isnan(dst[i+6]));


        const Real rInv = static_cast<Real>(1.0) / (r1 + r2);
        const Real alpha2 = max(static_cast<Real>(0.0), min(A2, static_cast<Real>(1.0)));
        const Real GmixInv = static_cast<Real>(1.0) / (g1m1Inv - alpha2*(g1m1Inv - g2m1Inv));
        Real Pmix = g1m1Inv*g1*pc1 - alpha2*(g1m1Inv*g1*pc1 - g2m1Inv*g2*pc2);

        const Real ke = static_cast<Real>(0.5)*rInv*(u*u + v*v + w*w);//whatever ke we had
        Real pressure = (e - ke - Pmix) * GmixInv;

        dst[i] = max(r1,static_cast<Real>(0.0));//change rho
        dst[i+1] = max(r2, static_cast<Real>(0.0));//change rho
        dst[i+6] = max(static_cast<Real>(0.0), min(A2, static_cast<Real>(1.0)));//change alpha2

#ifdef _NOK_
        if (Pmix*GmixInv/(GmixInv+1.0) < -m_alpha * pressure ) //if it was still bad with new P and new G
        {
            const Real difference = -m_beta * pressure*(GmixInv+1.0)/GmixInv - Pmix;
            Pmix += abs(difference);//change P again

            dst[i+5] = pressure/GmixInv + Pmix + ke;//update e
        }
#else
        // TODO: Ursula: I am still trying to figure out the best state clipping for
        // the volume-fraction-based 5eq system
        // January 2016: the simplest one currently performs best

        /*
           Real maxroot = (-g2*pc2 - alpha2 * (g1*pc1 - g2*pc2)) / (g2+alpha2*g1-alpha2*g2);
           const Real a = (pc1>0.0) ? -pc1 : 1.0e-6;
           const Real b = (pc2>0.0) ? -pc2 : 1.0e-6;
           if (alpha2 <= static_cast<Real>(0.0))
           maxroot = a;
           else if (alpha2 >= static_cast<Real>(1.0))
           maxroot = b;
           else
           maxroot =  max(static_cast<Real>(1.0e-6),max(maxroot,max(a,b)));
        //cout << "maxroot " << maxroot << " " << m_alpha  <<endl;

        assert(maxroot<=static_cast<Real>(1.0e-6));

        if (pressure <= maxroot/m_alpha)
        {
        //cout << "ENTER CORRECTION " << pressure << " " << pc1 << " " << pc2  << " " << maxroot << " " << alpha2 << endl;
        const Real difference = - m_beta * pressure + maxroot;
        pressure += abs(difference);
        // pressure += 1.0e-6;
        dst[i+5] = pressure/GmixInv + Pmix + ke;
        }
        */

        if (pressure<static_cast<Real>(1.0e-3))
        {
            pressure=1.0e-3;
            dst[i+5] = pressure/GmixInv + Pmix + ke;

            //cout << "CORRECTED PRESSURE" << endl;

        }

#ifndef NDEBUG
        const Real r1c1sq = g1*(pressure + pc1);
        const Real r2c2sq = g2*(pressure + pc2);
        const Real impedance = static_cast<Real>(1.0)/r1c1sq - alpha2 * (static_cast<Real>(1.0)/r1c1sq - static_cast<Real>(1.0)/r2c2sq);
        // if (impedance <= 0.0) cout << " IMPEDANCE " << impedance << " " << pressure << endl;
        assert(impedance>0.0);
#endif

        assert(dst[i+5] > 0);
#endif
    }
}
