/*
 *  Diffusion_CPP.cpp
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 2/27/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

#include "Diffusion_CPP.h"


// former, but in this form not necessary (see also header file)
void Diffusion_CPP::_convert(
    const Real * const gptfirst, const int gptfloats, const int rowgpts)
{
  InputSOA_ST &a=sa.ref();

  for(int sy=0; sy<_BLOCKSIZE_+2; sy++) { // raw indices (halo at 0)
    for(int sx=0; sx<_BLOCKSIZE_+2; sx++) {
      AssumedType pt = 
        *(AssumedType*)(gptfirst + gptfloats*(sx + sy*rowgpts));

      const int dx = sx-1;  // relative indices (halo at -1)
      const int dy = sy-1;

      a.ref(dx, dy) = pt.A2;
    }
  }
}


void Diffusion_CPP::_copyback(
    Real * const gptfirst, const int gptfloats, const int rowgpts)
{
  const Real c = dtinvh;
  for(int iy=0; iy<OutputSOA::NY; iy++) {
    for(int ix=0; ix<OutputSOA::NX; ix++) {
      AssumedType& rhs = 
        *(AssumedType*)(gptfirst + gptfloats*(ix + iy*rowgpts));
      rhs.A2  = rhsa(ix, iy);

      assert(isfinite(rhsa(ix, iy)));
    }
  }
}
