/*
 *  Diffusion_CPP.cpp
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 2/27/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

#include "Diffusion_CPP.h"


// b -- base pointer
// fe -- floats per element
// er -- elements per row
void Diffusion_CPP::_convert(const Real * const b, const int fe, const int er)
{
  InputSOA_ST &a=sa.ref();

  for(int sy=0; sy<_BLOCKSIZE_+2; sy++) { // abs indices (halo at 0)
    for(int sx=0; sx<_BLOCKSIZE_+2; sx++) {
      AssumedType pt = *(AssumedType*)(b + fe * (sx + sy * er));

      const int dx = sx-1;  // rel indices (halo at -1)
      const int dy = sy-1;

      a.ref(dx, dy) = pt.A2;
    }
  }
}


// b -- base pointer
// fe -- floats per element
// er -- elements per row
void Diffusion_CPP::_copyback(Real * const b, const int fe, const int er)
{
  const Real c = 0.1;
  for(int iy=0; iy<OutputSOA::NY; iy++) {
    for(int ix=0; ix<OutputSOA::NX; ix++) {
      AssumedType& rhs = *(AssumedType*)(b + fe * (ix + iy * er));
      rhs.A2  = rhsa(ix, iy) * c;

      assert(isfinite(rhsa(ix, iy)));
    }
  }
}
