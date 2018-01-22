/*
 *  DivTensor_CPP.cpp
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 2/24/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

#include <cassert>
#include <cmath>

#include "common.h"

#include "DivTensor_CPP.h"

//implemented as "conversion to primitive"
void DivTensor_CPP::_convert(const Real * const gptfirst, const int gptfloats, const int rowgpts)
{
    InputSOA_ST& u=ringu.ref(), &v=ringv.ref(), &w=ringw.ref(), &l = ringls.ref();

	for(int sy=0; sy<_BLOCKSIZE_+2; sy++)
		for(int sx=0; sx<_BLOCKSIZE_+2; sx++)
		{
            AssumedType pt = *(AssumedType*)(gptfirst + gptfloats*(sx + sy*rowgpts));

            const int dx = sx-1;
            const int dy = sy-1;

            const Real inv_r = static_cast<Real>(1.0)/(pt.r1 + pt.r2);

            u.ref(dx, dy) = inv_r * pt.u;
            v.ref(dx, dy) = inv_r * pt.v;
            w.ref(dx, dy) = inv_r * pt.w;

#if defined(_CONVERTCLIP_) || defined(_ALPHACLIP_)
            const Real alpha2 = min(static_cast<Real>(1.0), max(static_cast<Real>(0.0), pt.A2));
#else
            const Real alpha2 = pt.A2;
#endif /* defined(_CONVERTCLIP_) || defined(_ALPHACLIP_) */

            l.ref(dx, dy) = alpha2;
		}
}

void DivTensor_CPP::_corners(const InputSOA_ST& ls0, const InputSOA_ST& ls1, TempSOA_ST& nx, TempSOA_ST& ny, TempSOA_ST& nz)
{
	for(int iy=0; iy<TempSOA_ST::NY; iy++)
		for(int ix=0; ix<TempSOA_ST::NX; ix++)
		{
			const Real c000    = ls1(ix,iy);
			const Real cm100   = ls1(ix-1,iy);
			const Real c0m10   = ls1(ix,iy-1);
			const Real cm1m10  = ls1(ix-1,iy-1);
			const Real c00m1   = ls0(ix,iy);
			const Real cm10m1  = ls0(ix-1,iy);
			const Real c0m1m1  = ls0(ix,iy-1);
			const Real cm1m1m1 = ls0(ix-1,iy-1);

			nx.ref(ix, iy) = c000-cm100 + c0m10-cm1m10 + c00m1-cm10m1 + c0m1m1-cm1m1m1;
			ny.ref(ix, iy) = cm100-cm1m10 + c000-c0m10 + cm10m1-cm1m1m1  + c00m1-c0m1m1;
			nz.ref(ix, iy) = cm1m10-cm1m1m1 + cm100-cm10m1 + c0m10-c0m1m1  + c000-c00m1;
		}
}

//base tensor implementation is I*|n|/3 - nn^t
// Hint: for factor 1/3 see Hu & Adams 2006
void DivTensor_CPP::_tensor_xface(const TempSOA_ST& nx0, const TempSOA_ST& nx1,
								  const TempSOA_ST& ny0, const TempSOA_ST& ny1,
								  const TempSOA_ST& nz0, const TempSOA_ST& nz1)
{
#ifdef _SURFF3D_
	const Real factor = (Real)(1./3);
#elif _SURFF2D_
        const Real factor = (Real)(1./2);
#else
        const Real factor = (Real)(1.0);
#endif

	for(int iy=0; iy<TempPiXSOA_ST::NY; iy++)
        for(int ix=0; ix<TempPiXSOA_ST::NX; ix++)
        {
            const Real nx = (nx0(ix,iy)+nx0(ix,iy+1)+nx1(ix,iy)+nx1(ix,iy+1));
            const Real ny = (ny0(ix,iy)+ny0(ix,iy+1)+ny1(ix,iy)+ny1(ix,iy+1));
            const Real nz = (nz0(ix,iy)+nz0(ix,iy+1)+nz1(ix,iy)+nz1(ix,iy+1));

			const Real mag = mysqrt(nx*nx + ny*ny + nz*nz);
			const Real inv_mag = mag > static_cast<Real>(1.0e-9) ? ((Real)1)/mag : 0;
            assert(!isnan(inv_mag));

            txx.ref(ix, iy) = factor*mag - nx*nx*inv_mag;
			txy.ref(ix, iy) = -nx*ny*inv_mag;
			txz.ref(ix, iy) = -nx*nz*inv_mag;
		}
}

//base tensor implementation is I*|n|/3 - nn^t
void DivTensor_CPP::_tensor_yface(const TempSOA_ST& nx0, const TempSOA_ST& nx1,
								  const TempSOA_ST& ny0, const TempSOA_ST& ny1,
								  const TempSOA_ST& nz0, const TempSOA_ST& nz1)
{
#ifdef _SURFF3D_
	const Real factor = (Real)(1./3);
#elif _SURFF2D_
        const Real factor = (Real)(1./2);
#else
        const Real factor = (Real)(1.0);
#endif

	for(int iy=0; iy<TempPiYSOA_ST::NY; iy++)
        for(int ix=0; ix<TempPiYSOA_ST::NX; ix++)
        {
            const Real nx = (nx0(ix,iy)+nx0(ix+1,iy)+nx1(ix,iy)+nx1(ix+1,iy));
            const Real ny = (ny0(ix,iy)+ny0(ix+1,iy)+ny1(ix,iy)+ny1(ix+1,iy));
            const Real nz = (nz0(ix,iy)+nz0(ix+1,iy)+nz1(ix,iy)+nz1(ix+1,iy));

			const Real mag = mysqrt(nx*nx + ny*ny + nz*nz);
			const Real inv_mag = mag > static_cast<Real>(1.0e-9) ? ((Real)1)/mag : 0;
	    assert(!isnan(inv_mag));

            tyx.ref(ix, iy) = -ny*nx*inv_mag ;
			tyy.ref(ix, iy) = factor*mag - ny*ny*inv_mag;
			tyz.ref(ix, iy) = -ny*nz*inv_mag;
		}
}

//base tensor implementation is I*|n|/3 - nn^t
void DivTensor_CPP::_tensor_zface(const TempSOA_ST& nx0, const TempSOA_ST& ny0, const TempSOA_ST& nz0,
								  TempPiZSOA_ST& tzx, TempPiZSOA_ST& tzy, TempPiZSOA_ST& tzz)
{
#ifdef _SURFF3D_
	const Real factor = (Real)(1./3);
#elif _SURFF2D_
        const Real factor = (Real)(1./2);
#else
        const Real factor = (Real)(1.0);
#endif

	for(int iy=0; iy<TempPiZSOA_ST::NY; iy++)
        for(int ix=0; ix<TempPiZSOA_ST::NX; ix++)
        {
            const Real nx = (nx0(ix,iy)+nx0(ix+1,iy)+nx0(ix,iy+1)+nx0(ix+1,iy+1));
            const Real ny = (ny0(ix,iy)+ny0(ix+1,iy)+ny0(ix,iy+1)+ny0(ix+1,iy+1));
            const Real nz = (nz0(ix,iy)+nz0(ix+1,iy)+nz0(ix,iy+1)+nz0(ix+1,iy+1));

			const Real mag = mysqrt(nx*nx + ny*ny + nz*nz);
			const Real inv_mag = mag > static_cast<Real>(1.0e-9) ? ((Real)1)/mag : 0;
			assert(!isnan(inv_mag));

            tzx.ref(ix, iy) = -nz*nx*inv_mag ;
			tzy.ref(ix, iy) = -nz*ny*inv_mag;
			tzz.ref(ix, iy) = factor*mag - nz*nz*inv_mag;
		}
}

void DivTensor_CPP::_udot_tx(const InputSOA_ST& u, const InputSOA_ST& v, const InputSOA_ST& w)
{
	for(int iy=0; iy<TempPiXSOA_ST::NY; iy++)
		for(int ix=0; ix<TempPiXSOA_ST::NX; ix++)
		{
			utx.ref(ix, iy) =
			txx(ix, iy)*(u(ix, iy) + u(ix-1, iy)) +
			txy(ix, iy)*(v(ix, iy) + v(ix-1, iy)) +
			txz(ix, iy)*(w(ix, iy) + w(ix-1, iy));
		}
}

void DivTensor_CPP::_udot_ty(const InputSOA_ST& u, const InputSOA_ST& v, const InputSOA_ST& w)
{
	for(int iy=0; iy<TempPiYSOA_ST::NY; iy++)
		for(int ix=0; ix<TempPiYSOA_ST::NX; ix++)
		{
			uty.ref(ix, iy) =
			tyx(ix, iy)*(u(ix, iy) + u(ix, iy-1)) +
			tyy(ix, iy)*(v(ix, iy) + v(ix, iy-1)) +
			tyz(ix, iy)*(w(ix, iy) + w(ix, iy-1));
		}
}

void DivTensor_CPP::_udot_tz(const InputSOA_ST& u0, const InputSOA_ST& v0, const InputSOA_ST& w0,
							 const InputSOA_ST& u1, const InputSOA_ST& v1, const InputSOA_ST& w1,
							 const TempPiZSOA_ST& tzx, const TempPiZSOA_ST& tzy, const TempPiZSOA_ST& tzz,
							 TempPiZSOA_ST& utz)
{
	for(int iy=0; iy<TempPiZSOA_ST::NY; iy++)
		for(int ix=0; ix<TempPiZSOA_ST::NX; ix++)
		{
			utz.ref(ix, iy) =
			tzx(ix, iy)*(u0(ix, iy) + u1(ix, iy)) +
			tzy(ix, iy)*(v0(ix, iy) + v1(ix, iy)) +
			tzz(ix, iy)*(w0(ix, iy) + w1(ix, iy));
		}
}

void DivTensor_CPP::_div_dxy()
{
  for(int iy=0; iy<OutputSOA::NY; iy++)
    for(int ix=0; ix<OutputSOA::NX; ix++)
    {
      //rhsu.ref(ix, iy) = txx(ix+1, iy) - txx(ix, iy) + tyx(ix, iy+1) - tyx(ix, iy);
      //rhsv.ref(ix, iy) = txy(ix+1, iy) - txy(ix, iy) + tyy(ix, iy+1) - tyy(ix, iy);
      //rhsw.ref(ix, iy) = txz(ix+1, iy) - txz(ix, iy) + tyz(ix, iy+1) - tyz(ix, iy);
      //rhss.ref(ix, iy) = utx(ix+1, iy) - utx(ix, iy) + uty(ix, iy+1) - uty(ix, iy);
      rhsu.ref(ix, iy) = 0.;
      rhsv.ref(ix, iy) = 0.;
      rhsw.ref(ix, iy) = 0.;
      rhss.ref(ix, iy) = 0.;
    }
}

void DivTensor_CPP::_div_dz(const TempPiZSOA_ST& tzx0, const TempPiZSOA_ST& tzy0, const TempPiZSOA_ST& tzz0, const TempPiZSOA_ST& utz0,
							const TempPiZSOA_ST& tzx1, const TempPiZSOA_ST& tzy1, const TempPiZSOA_ST& tzz1, const TempPiZSOA_ST& utz1)
{
  for(int iy=0; iy<OutputSOA::NY; iy++)
    for(int ix=0; ix<OutputSOA::NX; ix++)
    {
      //rhsu.ref(ix, iy) += tzx1(ix, iy) - tzx0(ix, iy);
      //rhsv.ref(ix, iy) += tzy1(ix, iy) - tzy0(ix, iy);
      //rhsw.ref(ix, iy) += tzz1(ix, iy) - tzz0(ix, iy);
      //rhss.ref(ix, iy) += utz1(ix, iy) - utz0(ix, iy);
      rhsu.ref(ix, iy) = 0.;
      rhsv.ref(ix, iy) = 0.;
      rhsw.ref(ix, iy) = 0.;
      rhss.ref(ix, iy) = 0.;
    }
}

void DivTensor_CPP::_copyback(Real * const gptfirst, const int gptfloats, const int rowgpts)
{
	const Real factor = sigma*dtinvh/(16*h);

	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix++)
		{
			AssumedType& rhs = *(AssumedType*)(gptfirst + gptfloats*(ix + iy*rowgpts));
			rhs.u  = a*rhs.u + factor*rhsu(ix, iy);
			rhs.v  = a*rhs.v + factor*rhsv(ix, iy);
			rhs.w  = a*rhs.w + factor*rhsw(ix, iy);
            rhs.E  = a*rhs.E + factor*rhss(ix, iy)*0.5;

			assert(!isnan(rhsu(ix, iy)));
			assert(!isnan(rhsv(ix, iy)));
			assert(!isnan(rhsw(ix, iy)));
			assert(!isnan(rhss(ix, iy)));
		}
}
