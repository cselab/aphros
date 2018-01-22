/*
 *  MaxSpeedOfSound_CPP_5eq.cpp
 *  MPCFcore
 *
 *  Created by Fabian Wermelinger 10/22/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "common.h"
#include "MaxSpeedOfSound_CPP_5eq.h"


Real MaxSpeedOfSound_CPP_5eq::compute(const Real * const src, const int gptfloats) const
{
    const int N=_BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*gptfloats;
    Real sos = 0;

    for(int i=0; i<N; i+=gptfloats)
    {
        assert(!isnan(src[i]));
        assert(!isnan(src[i+1]));
        assert(!isnan(src[i+2]));
        assert(!isnan(src[i+3]));
        assert(!isnan(src[i+4]));
        assert(!isnan(src[i+5]));
        assert(!isnan(src[i+6]));

        assert(!isinf(src[i]));
        assert(!isinf(src[i+1]));
        assert(!isinf(src[i+2]));
        assert(!isinf(src[i+3]));
        assert(!isinf(src[i+4]));
        assert(!isinf(src[i+5]));
        assert(!isinf(src[i+6]));

        const Real ru = src[i+2];
        const Real rv = src[i+3];
        const Real rw = src[i+4];
        const Real E  = src[i+5];

#ifdef _CONVERTCLIP_
        const Real a1r1 = max(static_cast<Real>(0.0),src[i]);
        const Real a2r2 = max(static_cast<Real>(0.0),src[i+1]);
        const Real alpha2 = max(static_cast<Real>(ALPHAEPS), min(src[i+6], static_cast<Real>(1.0-ALPHAEPS)));
#else
        const Real a1r1 = src[i];
        const Real a2r2 = src[i+1];
#ifdef _ALPHACLIP_
        const Real alpha2 = max(static_cast<Real>(ALPHAEPS), min(src[i+6], static_cast<Real>(1.0-ALPHAEPS)));
#else
        const Real alpha2 = src[i+6];
#endif /* _ALPHACLIP_ */
#endif /* _CONVERTCLIP */
        const Real alpha1 = static_cast<Real>(1.0) - alpha2;

        const Real rInv = static_cast<Real>(1.0) / (a1r1 + a2r2);
        const Real GmixInv = static_cast<Real>(1.0) / (alpha1*_g1m1Inv + alpha2*_g2m1Inv);

#ifdef _CONVERTCLIP_
        const Real Pmix = max(static_cast<Real>(0.0),alpha1*_g1m1Inv*_g1*_pc1 + alpha2*_g2m1Inv*_g2*_pc2);
#else
        const Real Pmix = alpha1*_g1m1Inv*_g1*_pc1 + alpha2*_g2m1Inv*_g2*_pc2;
#endif /* _CONVERTCLIP_ */

        Real p = GmixInv*E - GmixInv*static_cast<Real>(0.5)*rInv*(ru*ru + rv*rv + rw*rw) - GmixInv*Pmix;

#ifdef _CONVERTCLIP_
#ifdef _NOK_
        const Real pthresh = Pmix*GmixInv/(GmixInv+1.0);
#else
        // TODO: (fabianw@mavt.ethz.ch; Mon 02 May 2016 02:46:04 PM CEST) What
        // to do here if ALPHAEPS is used?
        const Real pthresh = alpha2 == static_cast<Real>(1.0) ? _pc2 : (alpha2 == static_cast<Real>(0.0) ? _pc1 : min(_pc1, _pc2));
#endif
        const Real deltap = pthresh <= (-p + static_cast<Real>(PRESEPS)) ? (-p - pthresh + static_cast<Real>(PRESEPS)) : static_cast<Real>(0.0);
#ifdef _CLIPVERBOSE_
        if (deltap > static_cast<Real>(0.0))
        {
          cout << "Pressure correction in HLLC soskernel " << setprecision(12) << p << " " << deltap << " " << pthresh << " " << alpha1 << " " << alpha2 << endl;
        }
#endif
        p += deltap;
#endif /* _CONVERTCLIP_ */

#ifdef _NOK_
        const Real c = sqrt(((GmixInv + 1.0)*p + GmixInv*Pmix)*rInv);
#else
        // TODO: (fabianw@mavt.ethz.ch; Mon 02 May 2016 09:10:20 AM CEST)
        // testing
        const Real denom1 = _g1*(p + _pc1);
        const Real denom2 = _g2*(p + _pc2);
        Real rc2 = (denom1 > 0) ? alpha1/denom1 : 0;
        rc2 += (denom2 > 0) ? alpha2/denom2 : 0;
        const Real c2 = rInv/rc2;
        const Real c = sqrt(c2);

        // const Real c = sqrt(rInv*_g1*_g2*(p+_pc1)*(p+_pc2) / (alpha1*_g2*(p+_pc2) + alpha2*_g1*(p+_pc1)));
#endif

        assert(!isnan(p));
        assert(!isnan(c));
        assert(c > 0);

        sos = max(sos, c + max(max(abs(ru), abs(rv)), abs(rw))*rInv);
    }

    return sos;
}
