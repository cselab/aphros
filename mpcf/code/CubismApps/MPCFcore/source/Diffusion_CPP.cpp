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
// n -- floats per element
// row -- elements per row
void Diffusion_CPP::_convert(const Real * const b, const int n, const int row)
{
  InputSOA_ST &a=sa.ref();

  for(int sy=0; sy<_BLOCKSIZE_+2; sy++) { // raw indices (halo at 0)
    for(int sx=0; sx<_BLOCKSIZE_+2; sx++) {
      AssumedType pt = *(AssumedType*)(b + n * (sx + sy * row));

      const int dx = sx-1;  // relative indices (halo at -1)
      const int dy = sy-1;

      a.ref(dx, dy) = pt.A2;
    }
  }
}


// b -- base pointer
// n -- floats per element
// row -- elements per row
void Diffusion_CPP::_copyback(Real * const b, const int n, const int row)
{
  const Real c = dtinvh;
  for(int iy=0; iy<OutputSOA::NY; iy++) {
    for(int ix=0; ix<OutputSOA::NX; ix++) {
      AssumedType& rhs = *(AssumedType*)(b + n * (ix + iy * row));
      rhs.A2  = rhsa(ix, iy) * c;

      assert(isfinite(rhsa(ix, iy)));
    }
  }
}
