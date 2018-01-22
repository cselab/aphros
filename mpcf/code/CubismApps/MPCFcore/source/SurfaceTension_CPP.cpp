//
//  SurfaceTension_CPP.cpp
//  MPCFcore
//
//  Created by Babak Hejazialhosseini on 8/2/11.
//  Copyright 2011 ETH Zurich. All rights reserved.
//
#include <math.h>
#include <string.h>

#include "SurfaceTension_CPP.h"

void SurfaceTension_CPP::_convert(const Real * const gptfirst, const int gptfloats, const int rowgpts)
{
  InputSOA_ST& u=ringu.ref(), &v=ringv.ref(), &w=ringw.ref(), &l = ringls.ref();

  for(int sy=0; sy<_BLOCKSIZE_+2; sy++)
    for(int sx=0; sx<_BLOCKSIZE_+2; sx++)
      {
            AssumedType pt = *(AssumedType*)(gptfirst + gptfloats*(sx + sy*rowgpts));

            const int dx = sx-1;
            const int dy = sy-1;

#ifdef _CONVERTCLIP_
            const Real a1r1 = max(static_cast<Real>(0.0),pt.a1r1);
            const Real a2r2 = max(static_cast<Real>(0.0),pt.a2r2);
            const Real alpha2 = max(static_cast<Real>(ALPHAEPS), min(pt.A2, static_cast<Real>(1.0-ALPHAEPS)));
#else
            const Real a1r1 = pt.a1r1;
            const Real a2r2 = pt.a2r2;
#ifdef _ALPHACLIP_
            const Real alpha2 = max(static_cast<Real>(ALPHAEPS), min(pt.A2, static_cast<Real>(1.0-ALPHAEPS)));
#else
            const Real alpha2 = pt.A2;
#endif
#endif
            const Real rInv = static_cast<Real>(1.0) / (a1r1 + a2r2);

            u.ref(dx, dy) = pt.ru*rInv;
            v.ref(dx, dy) = pt.rv*rInv;
            w.ref(dx, dy) = pt.rw*rInv;

            // get liquid void fraction (see Perigaud and Saurel JCP 2005)
            l.ref(dx, dy) = alpha2;

            assert(!isnan(u.ref(dx, dy)));
            assert(!isnan(v.ref(dx, dy)));
            assert(!isnan(w.ref(dx, dy)));
            assert(!isnan(l.ref(dx, dy)));
        }
}
