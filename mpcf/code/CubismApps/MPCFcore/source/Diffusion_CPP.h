/*
 *  Diffusion_CPP.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 2/27/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */
#pragma once
#include <vector>
#include <cstdio>
#include <cmath>

#include "SOA2D.h"

typedef SOA2D<-1, _BLOCKSIZE_+1, -1, _BLOCKSIZE_+1> InputSOA_ST;
typedef RingSOA2D<-1, _BLOCKSIZE_+1, -1, _BLOCKSIZE_+1, 2> RingInputSOA_ST;

typedef SOA2D<0,_BLOCKSIZE_+1, 0,_BLOCKSIZE_+1> TempSOA_ST;
typedef RingSOA2D<0, _BLOCKSIZE_+1, 0,_BLOCKSIZE_+1, 2> RingTempSOA_ST;

class Diffusion_CPP
{
  protected:

    struct AssumedType { Real r1, r2, u, v, w, E, A2, dummy; };

    typedef Real RealTemp;

    void _convert(const Real * const gptfirst, const int gptfloats, const int rowgpts);
    void _copyback(Real* gptfirst, const int gptfloats, const int rowgpts);
    Real dtinvh;

  public:
    RingInputSOA_ST sa; // input slice (source)
    OutputSOA rhsa; // output slice (rhs)

    Diffusion_CPP(const Real a,
        const Real mu1,
        const Real mu2,
        const Real h,
        const Real smoothing_length,
        const Real dtinvh) 
      : dtinvh(dtinvh) 
    {} 

    void compute(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
        Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts)
    {
      _convert(srcfirst, srcfloats, rowsrcs);
      sa.next();
      _convert(srcfirst + srcfloats*slicesrcs, srcfloats, rowsrcs);

      for(int islice=0; islice<_BLOCKSIZE_; islice++)
      {
        sa.next();

        _convert(srcfirst + (islice+2)*srcfloats*slicesrcs, srcfloats, rowsrcs);

        for(int iy=0; iy<OutputSOA::NY; ++iy) {
          for(int ix=0; ix<OutputSOA::NX; ++ix) {
            //rhsa.ref(ix, iy) = txx(ix+1, iy) - txx(ix, iy) + tyx(ix, iy+1) - tyx(ix, iy);
            rhsa.ref(ix, iy) = 0.;
          }
        }

        _copyback(dstfirst + islice*dstfloats*slicedsts, dstfloats, rowdsts);
      }
    }

    void hpc_info(float& flop_convert, int& traffic_convert,
        float& flop_corners, int& traffic_corner,
        float& flop_tensorface, int& traffic_tensorface,
        float& flop_tdotu, int& traffic_tdotu,
        float& flop_div, int& traffic_div,
        float& flop_copyback, int& traffic_copyback,
        size_t& footprint)
    {

      const int ninputs = (int)powf(_BLOCKSIZE_ + 2, 3);
      const int ntensors = 3 * (_BLOCKSIZE_ + 1) * _BLOCKSIZE_ * _BLOCKSIZE_;

      flop_convert = 13 * ninputs;
      traffic_convert = (4 + 4) * sizeof(Real) * ninputs;

      flop_corners *= 3;
      traffic_corner *= 3;

      flop_tensorface =  32 * ntensors;
      traffic_tensorface = (30 + 3) * sizeof(Real) * ntensors;

      footprint = sizeof(*this);
    }
};
