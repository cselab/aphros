/*
 *  Diffusion_CPP.cpp
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 2/27/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

#include "Diffusion_CPP.h"

// unnecessary (see also header file)
/*
void Diffusion_CPP::_convert_single_phase(const Real * const gptfirst, const int gptfloats, const int rowgpts)
{
    InputSOA_ST &u=ringu.ref(), &v=ringv.ref(), &w=ringw.ref(), &mu=ringls.ref();

    for(int sy=0; sy<_BLOCKSIZE_+2; sy++)
        for(int sx=0; sx<_BLOCKSIZE_+2; sx++)
        {
            AssumedType pt = *(AssumedType*)(gptfirst + gptfloats*(sx + sy*rowgpts));

            const int dx = sx-1;
            const int dy = sy-1;

            // This assumes that the phase is described by the index "1"
            // variables
            const Real inv_r = static_cast<Real>(1.0)/pt.r1;

            u.ref(dx, dy) = inv_r * pt.u;
            v.ref(dx, dy) = inv_r * pt.v;
            w.ref(dx, dy) = inv_r * pt.w;
            mu.ref(dx, dy) = mu1;
        }
}
*/

// former, but in this form not necessary (see also header file)
//void Diffusion_CPP::_convert_two_phase(const Real * const gptfirst, const int gptfloats, const int rowgpts)
void Diffusion_CPP::_convert(const Real * const gptfirst, const int gptfloats, const int rowgpts)
{
    InputSOA_ST &u=ringu.ref(), &v=ringv.ref(), &w=ringw.ref(), &mu=ringls.ref();

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
            const Real alpha1 = static_cast<Real>(1.0) - alpha2;

            mu.ref(dx, dy) = alpha1*mu1 + alpha2*mu2;
        }
}


void Diffusion_CPP::_xface(const TempSOA_ST& ux0, const TempSOA_ST& uy0, const TempSOA_ST& uz0, const TempSOA_ST& ux1, const TempSOA_ST& uy1, const TempSOA_ST& uz1,
        const TempSOA_ST& vx0, const TempSOA_ST& vy0, const TempSOA_ST& vx1, const TempSOA_ST& vy1,
        const TempSOA_ST& wx0, const TempSOA_ST& wz0, const TempSOA_ST& wx1, const TempSOA_ST& wz1,
        const InputSOA_ST &mu0)
{
    const RealTemp factor = static_cast<RealTemp>(2.0/3.0);

    for(int iy=0; iy<TempPiXSOA_ST::NY; iy++)
        for(int ix=0; ix<TempPiXSOA_ST::NX; ix++)
        {
            const RealTemp dudx = (ux0(ix,iy)+ux0(ix,iy+1)+ux1(ix,iy)+ux1(ix,iy+1));
            const RealTemp dudy = (uy0(ix,iy)+uy0(ix,iy+1)+uy1(ix,iy)+uy1(ix,iy+1));
            const RealTemp dudz = (uz0(ix,iy)+uz0(ix,iy+1)+uz1(ix,iy)+uz1(ix,iy+1));
            const RealTemp dvdx = (vx0(ix,iy)+vx0(ix,iy+1)+vx1(ix,iy)+vx1(ix,iy+1));
            const RealTemp dvdy = (vy0(ix,iy)+vy0(ix,iy+1)+vy1(ix,iy)+vy1(ix,iy+1));
            const RealTemp dwdx = (wx0(ix,iy)+wx0(ix,iy+1)+wx1(ix,iy)+wx1(ix,iy+1));
            const RealTemp dwdz = (wz0(ix,iy)+wz0(ix,iy+1)+wz1(ix,iy)+wz1(ix,iy+1));

            const RealTemp mu = static_cast<RealTemp>(0.5) * (mu0(ix-1,iy) + mu0(ix,iy));

            txx.ref(ix, iy) = mu * (static_cast<RealTemp>(2.0)*dudx - factor * (dudx + dvdy + dwdz));
            txy.ref(ix, iy) = mu * (dvdx + dudy);
            txz.ref(ix, iy) = mu * (dwdx + dudz);
        }
}

void Diffusion_CPP::_yface(const TempSOA_ST& ux0, const TempSOA_ST& uy0, const TempSOA_ST& ux1, const TempSOA_ST& uy1,
        const TempSOA_ST& vx0, const TempSOA_ST& vy0, const TempSOA_ST& vz0, const TempSOA_ST& vx1, const TempSOA_ST& vy1, const TempSOA_ST& vz1,
        const TempSOA_ST& wy0, const TempSOA_ST& wz0, const TempSOA_ST& wy1, const TempSOA_ST& wz1,
        const InputSOA_ST &mu0)
{
    const RealTemp factor = static_cast<RealTemp>(2.0/3.0);

    for(int iy=0; iy<TempPiYSOA_ST::NY; iy++)
        for(int ix=0; ix<TempPiYSOA_ST::NX; ix++)
        {
            const RealTemp dudx = (ux0(ix,iy)+ux0(ix+1,iy)+ux1(ix,iy)+ux1(ix+1,iy));
            const RealTemp dudy = (uy0(ix,iy)+uy0(ix+1,iy)+uy1(ix,iy)+uy1(ix+1,iy));
            const RealTemp dvdx = (vx0(ix,iy)+vx0(ix+1,iy)+vx1(ix,iy)+vx1(ix+1,iy));
            const RealTemp dvdy = (vy0(ix,iy)+vy0(ix+1,iy)+vy1(ix,iy)+vy1(ix+1,iy));
            const RealTemp dvdz = (vz0(ix,iy)+vz0(ix+1,iy)+vz1(ix,iy)+vz1(ix+1,iy));
            const RealTemp dwdy = (wy0(ix,iy)+wy0(ix+1,iy)+wy1(ix,iy)+wy1(ix+1,iy));
            const RealTemp dwdz = (wz0(ix,iy)+wz0(ix+1,iy)+wz1(ix,iy)+wz1(ix+1,iy));

            const RealTemp mu = static_cast<RealTemp>(0.5) * (mu0(ix,iy-1) + mu0(ix,iy));

            tyx.ref(ix, iy) = mu * (dudy + dvdx);
            tyy.ref(ix, iy) = mu * (static_cast<RealTemp>(2.0)*dvdy - factor*(dudx + dvdy + dwdz));
            tyz.ref(ix, iy) = mu * (dwdy + dvdz);
        }
}

void Diffusion_CPP::_zface(const TempSOA_ST& ux, const TempSOA_ST& uz,
        const TempSOA_ST& vy, const TempSOA_ST& vz,
        const TempSOA_ST& wx, const TempSOA_ST& wy, const TempSOA_ST& wz,
        const InputSOA_ST &mu0, const InputSOA_ST &mu1,
        TempPiZSOA_ST& tzx, TempPiZSOA_ST& tzy, TempPiZSOA_ST& tzz)
{
    const RealTemp factor = static_cast<RealTemp>(2.0/3.0);

    for(int iy=0; iy<TempPiZSOA_ST::NY; iy++)
        for(int ix=0; ix<TempPiZSOA_ST::NX; ix++)
        {
            const RealTemp dudx = (ux(ix,iy)+ux(ix+1,iy)+ux(ix,iy+1)+ux(ix+1,iy+1));
            const RealTemp dudz = (uz(ix,iy)+uz(ix+1,iy)+uz(ix,iy+1)+uz(ix+1,iy+1));
            const RealTemp dvdy = (vy(ix,iy)+vy(ix+1,iy)+vy(ix,iy+1)+vy(ix+1,iy+1));
            const RealTemp dvdz = (vz(ix,iy)+vz(ix+1,iy)+vz(ix,iy+1)+vz(ix+1,iy+1));
            const RealTemp dwdx = (wx(ix,iy)+wx(ix+1,iy)+wx(ix,iy+1)+wx(ix+1,iy+1));
            const RealTemp dwdy = (wy(ix,iy)+wy(ix+1,iy)+wy(ix,iy+1)+wy(ix+1,iy+1));
            const RealTemp dwdz = (wz(ix,iy)+wz(ix+1,iy)+wz(ix,iy+1)+wz(ix+1,iy+1));

            const RealTemp mu = static_cast<RealTemp>(0.5) * (mu0(ix,iy) + mu1(ix,iy));

            tzx.ref(ix, iy) = mu * (dudz + dwdx);
            tzy.ref(ix, iy) = mu * (dvdz + dwdy);
            tzz.ref(ix, iy) = mu * (static_cast<RealTemp>(2.0)*dwdz - factor*(dudx + dvdy + dwdz));
        }
}
