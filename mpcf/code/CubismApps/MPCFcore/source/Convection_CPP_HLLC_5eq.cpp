/*
 *  Convection_CPP_HLLC_5eq.cpp
 *  MPCFcore
 *
 *  Created by Fabian Wermelinger 10/21/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#include <cmath>
#include <cassert>
#include <iostream>
#include <cstdio>

using namespace std;

#include "Convection_CPP_HLLC_5eq.h"

// #define LAMBDA

void Convection_CPP_HLLC_5eq::_convert(const Real * const gptfirst, const int gptfloats, const int rowgpts, const int islice=0)
{
    InputSOA& rho1 = this->rho1.ring.ref();
    InputSOA& rho2 = this->rho2.ring.ref();
    InputSOA& rho = this->rho.ring.ref();
    InputSOA& u    = this->u.ring.ref();
    InputSOA& v    = this->v.ring.ref();
    InputSOA& w    = this->w.ring.ref();
    InputSOA& p    = this->p.ring.ref();
    InputSOA& A2   = this->A2.ring.ref();

    for(int sy=0; sy<_BLOCKSIZE_+6; sy++)
        for(int sx=0; sx<_BLOCKSIZE_+6; sx++)
        {
            AssumedType pt = *(AssumedType*)(gptfirst + gptfloats*(sx + sy*rowgpts));

            const int dx = sx-3;
            const int dy = sy-3;
/*
           assert(!std::isinf(pt.r1));   assert(!std::isnan(pt.r1));
           assert(!std::isinf(pt.r2));   assert(!std::isnan(pt.r2));
           assert(!std::isinf(pt.u));    assert(!std::isnan(pt.u));
           assert(!std::isinf(pt.v));    assert(!std::isnan(pt.v));
           assert(!std::isinf(pt.w));    assert(!std::isnan(pt.w));
           assert(!std::isinf(pt.E));    assert(!std::isnan(pt.E));
           assert(!std::isinf(pt.A2));   assert(!std::isnan(pt.A2));
*/

#ifdef _CONVERTCLIP_
            const Real a1r1 = max(static_cast<Real>(0.0),pt.r1);
            const Real a2r2 = max(static_cast<Real>(0.0),pt.r2);
            const Real alpha2 = max(static_cast<Real>(ALPHAEPS), min(pt.A2, static_cast<Real>(1.0-ALPHAEPS)));
            const Real alpha1 = static_cast<Real>(1.0)-alpha2;

            const Real rInv = static_cast<Real>(1.0) / (a1r1 + a2r2);
            const Real GmixInv = static_cast<Real>(1.0) / (alpha1*_g1m1Inv + alpha2*_g2m1Inv);
            const Real Pmix = max(static_cast<Real>(0.0),alpha1*_g1m1Inv*_g1*_pc1 + alpha2*_g2m1Inv*_g2*_pc2);

            /* TODO: (fabianw; Wed 21 Oct 2015 09:28:19 PM CEST) check relative
             * diff of GmixInv, Pmix and the result of the subtraction in
             * the parenthesis! (Ursula: also affects sos kernel) */
            /* p.ref(dx, dy) = GmixInv*(pt.E - static_cast<Real>(0.5)*rInv*(pt.u*pt.u + pt.v*pt.v + pt.w*pt.w) - Pmix); */
            const Real actpres = GmixInv*pt.E - GmixInv*static_cast<Real>(0.5)*rInv*(pt.u*pt.u + pt.v*pt.v + pt.w*pt.w) - GmixInv*Pmix;
#ifdef _NOK_
            const Real pthresh = Pmix*GmixInv/(GmixInv+1.0);
#else
            // TODO: (fabianw@mavt.ethz.ch; Mon 02 May 2016 02:35:33 PM CEST)
            // what to do here if ALPHAEPS is used?
            const Real pthresh = alpha2 == static_cast<Real>(1.0) ? _pc2 : (alpha2 == static_cast<Real>(0.0) ? _pc1 : min(_pc1, _pc2));
#endif
//            const Real deltap = pthresh < (-static_cast<Real>(2.0)*actpres) ? (-static_cast<Real>(4.0)*actpres - pthresh) : static_cast<Real>(0.0);
            const Real deltap =  pthresh <= (-actpres + static_cast<Real>(PRESEPS)) ? (-actpres - pthresh + static_cast<Real>(PRESEPS)) : static_cast<Real>(0.0);
/*
#ifdef _CLIPVERBOSE_
            if (deltap > static_cast<Real>(0.0))
            {
              cout << "Pressure correction in HLLC convert " << setprecision(12) << actpres << " " << deltap << " " << pthresh << " " << alpha1 << " " << alpha2  << endl;
            }
#endif
*/
            const Real pressure = actpres + deltap;
#ifdef _QUANTIFY_DELTAP_
            if (deltap > 0)
            {
                _pcorrMax = (deltap > _pcorrMax) ? deltap : _pcorrMax;
                ++_pcorrCount;
            }
#endif /* _QUANTIFY_DELTAP_ */

            rho1.ref(dx, dy) = a1r1;
            rho2.ref(dx, dy) = a2r2;
            u.ref(dx, dy) = pt.u*rInv;
            v.ref(dx, dy) = pt.v*rInv;
            w.ref(dx, dy) = pt.w*rInv;
            p.ref(dx, dy) = pressure;
            A2.ref(dx, dy) = alpha2;
#ifdef _A2_A2R2_R_
            rho.ref(dx, dy) = a1r1 + a2r2;
#endif

#else

            const Real rInv = static_cast<Real>(1.0) / (pt.r1 + pt.r2);
#ifdef _ALPHACLIP_
            const Real alpha2 = max(static_cast<Real>(ALPHAEPS), min(pt.A2, static_cast<Real>(1.0-ALPHAEPS)));
#else
            const Real alpha2 = pt.A2;
#endif /* _ALPHACLIP_ */
            const Real alpha1 = static_cast<Real>(1.0) - alpha2;
            const Real GmixInv = static_cast<Real>(1.0) / (alpha1*_g1m1Inv + alpha2*_g2m1Inv);
            const Real Pmix  = alpha1*_g1m1Inv*_g1*_pc1 + alpha2*_g2m1Inv*_g2*_pc2;
            const Real pres = GmixInv*pt.E - GmixInv*static_cast<Real>(0.5)*rInv*(pt.u*pt.u + pt.v*pt.v + pt.w*pt.w) - GmixInv*Pmix;
            rho1.ref(dx, dy) = pt.r1;
            rho2.ref(dx, dy) = pt.r2;
            u.ref(dx, dy) = pt.u*rInv;
            v.ref(dx, dy) = pt.v*rInv;
            w.ref(dx, dy) = pt.w*rInv;
            p.ref(dx, dy) = pres;
            A2.ref(dx, dy) = pt.A2;
#ifdef _A2_A2R2_R_
            rho.ref(dx, dy) = pt.r1+pt.r2;
#endif

#endif

/* #ifndef NDEBUG */
/*             const bool bEdgeXY = (sx<3 || sx>(_BLOCKSIZE_+2)) && (sy<3 || sy>(_BLOCKSIZE_+2)); */
/*             const bool bEdgeYZ = (sy<3|| sy>(_BLOCKSIZE_+2)) && (islice<3 || islice>(_BLOCKSIZE_+2)); */
/*             const bool bEdgeZX = (islice<3 || islice>(_BLOCKSIZE_+2)) && (sx<3  || sx>(_BLOCKSIZE_+2)); */
/*             const bool bCorner = (sx<3 || sx>(_BLOCKSIZE_+2)) && (sy<3 || sy>(_BLOCKSIZE_+2)) && (islice<3 || islice>(_BLOCKSIZE_+2)); */

/*             const bool bTensorialZone = bEdgeXY || bEdgeYZ || bEdgeZX || bCorner; */

/*             if (!bTensorialZone) */
/*             { */
/*                 if(std::isnan(pt.r1)) */
/*                     printf("nan pt.r1 at %d %d\n", dx, dy); */

/*                 if(std::isnan(pt.r2)) */
/*                     printf("nan pt.r2 at %d %d\n", dx, dy); */

/*                 if(pt.r1<=0.0) */
/*                     printf("negative pt.r1 at %d %d %e\n", dx, dy, pt.r1); */

/*                 if(pt.r2<=0.0) */
/*                     printf("negative pt.r2 at %d %d %e\n", dx, dy, pt.r2); */

/*                 if(pt.A2<0.0) */
/*                     printf("negative pt.A2 at %d %d %e\n", dx, dy, pt.A2); */
/*             } */
/* #endif */
        }
}

#ifdef _RECONPCLIP_
void Convection_CPP_HLLC_5eq::_pressure_clipping_hllc(const int SizeX,
        const TempSOA& A2m, const TempSOA& A2p,
        TempSOA& pm, TempSOA& pp)
{
    for (int iy = 0; iy < TempSOA::NY; ++iy)
        for (int ix = 0; ix < SizeX; ++ix)
        {
            // NOTE: This function is called if _RECONPCLIP_ is set at compile
            // time and does always clip, no matter if _CONVERTCLIP_ or
            // _ALPHACLIP_ are set.
            const Real alpha2m = max(static_cast<Real>(ALPHAEPS), min(A2m(ix,iy), static_cast<Real>(1.0-ALPHAEPS)));
            const Real alpha2p = max(static_cast<Real>(ALPHAEPS), min(A2p(ix,iy), static_cast<Real>(1.0-ALPHAEPS)));
            const Real alpha1m = 1.0 - alpha2m;
            const Real alpha1p = 1.0 - alpha2p;

#ifdef _NOK_
            const Real GmixInvm = static_cast<Real>(1.0) / (alpha1m*_g1m1Inv + alpha2m*_g2m1Inv);
            const Real GmixInvp = static_cast<Real>(1.0) / (alpha1p*_g1m1Inv + alpha2p*_g2m1Inv);
            const Real Pmixm = max(static_cast<Real>(0.0), alpha1m*_g1m1Inv*_g1*_pc1 + alpha2m*_g2m1Inv*_g2*_pc2);
            const Real Pmixp = max(static_cast<Real>(0.0), alpha1p*_g1m1Inv*_g1*_pc1 + alpha2p*_g2m1Inv*_g2*_pc2);

            const Real pthreshm = Pmixm*GmixInvm/(GmixInvm+1.0);
            const Real pthreshp = Pmixp*GmixInvp/(GmixInvp+1.0);

#else
            const Real pthreshm = alpha2m == static_cast<Real>(1.0) ? _pc2 : (alpha2m == static_cast<Real>(0.0) ? _pc1 : min(_pc1, _pc2));
            const Real pthreshp = alpha2p == static_cast<Real>(1.0) ? _pc2 : (alpha2p == static_cast<Real>(0.0) ? _pc1 : min(_pc1, _pc2));
#endif

//            const Real deltapm = pthreshm < (-static_cast<Real>(2.0)*pm(ix,iy)) ? (-static_cast<Real>(4.0)*pm(ix,iy) - pthreshm) : static_cast<Real>(0.0);
//            const Real deltapp = pthreshp < (-static_cast<Real>(2.0)*pp(ix,iy)) ? (-static_cast<Real>(4.0)*pp(ix,iy) - pthreshp) : static_cast<Real>(0.0);

            // Always clipped
            const Real deltapm =  pthreshm <= (-pm(ix,iy) + static_cast<Real>(PRESEPS)) ? (-pm(ix,iy) - pthreshm + static_cast<Real>(PRESEPS)) : static_cast<Real>(0.0);
            const Real deltapp =  pthreshp <= (-pp(ix,iy) + static_cast<Real>(PRESEPS)) ? (-pp(ix,iy) - pthreshp + static_cast<Real>(PRESEPS)) : static_cast<Real>(0.0);
#ifdef _CLIPVERBOSE_
            if (deltapm > static_cast<Real>(0.0))
            {
              cout << "Pressure correction in HLLC reconstruction " << setprecision(12) << pm(ix,iy) << " " << deltapm << " " << pthreshm << " " << alpha1m << " " << alpha2m  << endl;
            }
            if (deltapp > static_cast<Real>(0.0))
            {
              cout << "Pressure correction in HLLC reconstruction " << setprecision(12) << pp(ix,iy) << " " << deltapp << " " << pthreshp << " " << alpha1p << " " << alpha2p  << endl;
            }
#endif
            pm.ref(ix,iy) += deltapm;
            pp.ref(ix,iy) += deltapp;
        }
}
#endif /* _RECONPCLIP_ */

/* *
 * Right hand side computation for the transport equations of the advected
 * quantities.  The computation of the velocity at the cell interface follows
 * E. Johnsen and T. Colonius, "Implementation of WENO schemes in compressible
 * multicomponent flow problems.", Journal of Computational Physics, (2006).
 * */
void Convection_CPP_HLLC_5eq::_xextraterm_hllc( const TempSOA& um, const TempSOA& up,
        const TempSOA& A2m, const TempSOA& A2p,
        const TempSOA& pm, const TempSOA& pp,
        const TempSOA& sm, const TempSOA& sp, const TempSOA& ss)
{
    for (int iy = 0; iy < OutputSOA::NY; ++iy)
        for (int ix = 0; ix < OutputSOA::NX; ++ix)
        {
            int sign_star;
            Real s_minus, s_pluss, chi_starm, chi_starp;

            // right cell interface (u_{i+1/2})
            sign_star = ss(ix+1, iy) == 0 ? 0 : (ss(ix+1, iy) < 0 ? -1 : 1);
            s_minus   = min(static_cast<Real>(0.0), sm(ix+1, iy));
            s_pluss   = max(static_cast<Real>(0.0), sp(ix+1, iy));
            chi_starm = (sm(ix+1, iy) - um(ix+1, iy))/(sm(ix+1, iy) - ss(ix+1, iy)) - static_cast<Real>(1.0);
            chi_starp = (sp(ix+1, iy) - up(ix+1, iy))/(sp(ix+1, iy) - ss(ix+1, iy)) - static_cast<Real>(1.0);
            const Real u_hllc1 = static_cast<Real>(0.5*(1 + sign_star))*(um(ix+1, iy) + s_minus*chi_starm) + static_cast<Real>(0.5*(1 - sign_star))*(up(ix+1, iy) + s_pluss*chi_starp);

            // left cell interface (u_{i-1/2})
            sign_star = ss(ix, iy) == 0 ? 0 : (ss(ix, iy) < 0 ? -1 : 1);
            s_minus   = min(static_cast<Real>(0.0), sm(ix, iy));
            s_pluss   = max(static_cast<Real>(0.0), sp(ix, iy));
            chi_starm = (sm(ix, iy) - um(ix, iy))/(sm(ix, iy) - ss(ix, iy)) - static_cast<Real>(1.0);
            chi_starp = (sp(ix, iy) - up(ix, iy))/(sp(ix, iy) - ss(ix, iy)) - static_cast<Real>(1.0);
            const Real u_hllc0 = static_cast<Real>(0.5*(1 + sign_star))*(um(ix, iy) + s_minus*chi_starm) + static_cast<Real>(0.5*(1 - sign_star))*(up(ix, iy) + s_pluss*chi_starp);

            divu.ref(ix, iy) = u_hllc1 - u_hllc0;
        }

    for (int iy = 0; iy < OutputSOA::NY; ++iy)
        for (int ix = 0; ix < OutputSOA::NX; ++ix)
        {
            const Real a2p = A2p(ix,iy);
            const Real a2m = A2m(ix+1,iy);
            const Real a1p = static_cast<Real>(1.0) - a2p;
            const Real a1m = static_cast<Real>(1.0) - a2m;

            // A2
            //sumA2.ref(ix, iy) = A2p(ix, iy) + A2m(ix+1, iy);
            sumA2.ref(ix, iy) = a2p + a2m;

            // K
#ifndef LAMBDA
            // TODO: (fabianw@mavt.ethz.ch; Mon 02 May 2016 09:26:04 AM CEST)
            // testing
            const Real num1m = max(static_cast<Real>(0.0), _g1*(pm(ix+1,iy) + _pc1));
            const Real num2m = max(static_cast<Real>(0.0), _g2*(pm(ix+1,iy) + _pc2));
            const Real num1p = max(static_cast<Real>(0.0), _g1*(pp(ix,iy) + _pc1));
            const Real num2p = max(static_cast<Real>(0.0), _g2*(pp(ix,iy) + _pc2));
            sumK.ref(ix,iy) = 0;
            // sumK.ref(ix,iy) += (0 < a1m && a1m < 1) ? (num1m - num2m)/(num1m/a1m + num2m/a2m) : 0;
            // sumK.ref(ix,iy) += (0 < a1p && a1p < 1) ? (num1p - num2p)/(num1p/a1p + num2p/a2p) : 0;
            sumK.ref(ix,iy) += (static_cast<Real>(ALPHAEPS) < a1m && a1m < static_cast<Real>(1-ALPHAEPS)) ? (num1m - num2m)/(num1m/a1m + num2m/a2m) : 0;
            sumK.ref(ix,iy) += (static_cast<Real>(ALPHAEPS) < a1p && a1p < static_cast<Real>(1-ALPHAEPS)) ? (num1p - num2p)/(num1p/a1p + num2p/a2p) : 0;

            // sumK.ref(ix,iy) = a1p*a2p*(pp(ix,iy)*(_g1-_g2) + _g1*_pc1 - _g2*_pc2)/(pp(ix,iy)*(a1p*_g2 + a2p*_g1) + a1p*_g2*_pc2 + a2p*_g1*_pc1)
            // + a1m*a2m*(pm(ix+1,iy)*(_g1-_g2) + _g1*_pc1 - _g2*_pc2)/(pm(ix+1,iy)*(a1m*_g2 + a2m*_g1) + a1m*_g2*_pc2 + a2m*_g1*_pc1);
#else
            const Real gpc1m = _g1*(pm(ix+1,iy) + _pc1);
            const Real gpc2m = _g2*(pm(ix+1,iy) + _pc2);
            const Real gpc1p = _g1*(pp(ix,iy) + _pc1);
            const Real gpc2p = _g2*(pp(ix,iy) + _pc2);

            const Real sela2m = (a2m > 0) ? 1.0 : -1.0;
            const Real sela2p = (a2p > 0) ? 1.0 : -1.0;
            const Real selgpc1m = (gpc1m > 0) ? 1.0 : -1.0;
            const Real selgpc2m = (gpc2m > 0) ? 1.0 : -1.0;
            const Real selgpc1p = (gpc1p > 0) ? 1.0 : -1.0;
            const Real selgpc2p = (gpc2p > 0) ? 1.0 : -1.0;

            const Real a2m_pos   = (sela2m < 0)   ? 1.0 : a2m;
            const Real gpc1m_pos = (selgpc1m < 0) ? 1.0 : gpc1m;
            const Real gpc2m_pos = (selgpc2m < 0) ? 0.0 : gpc2m;
            const Real a2p_pos   = (sela2p < 0)   ? 1.0 : a2p;
            const Real gpc1p_pos = (selgpc1p < 0) ? 1.0 : gpc1p;
            const Real gpc2p_pos = (selgpc2p < 0) ? 0.0 : gpc2p;

            const Real a1m_pos = 1.0 - a2m_pos;
            const Real a1p_pos = 1.0 - a2p_pos;
            const Real lambdam = 0.25*(1.0+sela2m)*(1.0+selgpc1m);
            const Real lambdap = 0.25*(1.0+sela2p)*(1.0+selgpc1p);

            sumK.ref(ix,iy) = 0;
            sumK.ref(ix,iy) += lambdam*a1m_pos*a2m_pos*(gpc1m_pos - gpc2m_pos) / (a1m_pos*gpc2m_pos + a2m_pos*gpc1m_pos);
            sumK.ref(ix,iy) += lambdap*a1p_pos*a2p_pos*(gpc1p_pos - gpc2p_pos) / (a1p_pos*gpc2p_pos + a2p_pos*gpc1p_pos);
#endif /* 0 */
        }
}

void Convection_CPP_HLLC_5eq::_yextraterm_hllc( const TempSOA& vm, const TempSOA& vp,
        const TempSOA& A2m, const TempSOA& A2p,
        const TempSOA& pm, const TempSOA& pp,
        const TempSOA& sm, const TempSOA& sp, const TempSOA& ss)
{
    for (int iy = 0; iy < OutputSOA::NY; ++iy)
        for (int ix = 0; ix < OutputSOA::NX; ++ix)
        {
            int sign_star;
            Real s_minus, s_pluss, chi_starm, chi_starp;

            // right cell interface (v_{i+1/2})
            sign_star = ss(iy+1, ix) == 0 ? 0 : (ss(iy+1, ix) < 0 ? -1 : 1);
            s_minus   = min(static_cast<Real>(0.0), sm(iy+1, ix));
            s_pluss   = max(static_cast<Real>(0.0), sp(iy+1, ix));
            chi_starm = (sm(iy+1, ix) - vm(iy+1, ix))/(sm(iy+1, ix) - ss(iy+1, ix)) - static_cast<Real>(1.0);
            chi_starp = (sp(iy+1, ix) - vp(iy+1, ix))/(sp(iy+1, ix) - ss(iy+1, ix)) - static_cast<Real>(1.0);
            const Real v_hllc1 = static_cast<Real>(0.5*(1 + sign_star))*(vm(iy+1, ix) + s_minus*chi_starm) + static_cast<Real>(0.5*(1 - sign_star))*(vp(iy+1, ix) + s_pluss*chi_starp);

            // left cell interface (v_{i-1/2})
            sign_star = ss(iy, ix) == 0 ? 0 : (ss(iy, ix) < 0 ? -1 : 1);
            s_minus   = min(static_cast<Real>(0.0), sm(iy, ix));
            s_pluss   = max(static_cast<Real>(0.0), sp(iy, ix));
            chi_starm = (sm(iy, ix) - vm(iy, ix))/(sm(iy, ix) - ss(iy, ix)) - static_cast<Real>(1.0);
            chi_starp = (sp(iy, ix) - vp(iy, ix))/(sp(iy, ix) - ss(iy, ix)) - static_cast<Real>(1.0);
            const Real v_hllc0 = static_cast<Real>(0.5*(1 + sign_star))*(vm(iy, ix) + s_minus*chi_starm) + static_cast<Real>(0.5*(1 - sign_star))*(vp(iy, ix) + s_pluss*chi_starp);

            divu.ref(ix, iy) += v_hllc1 - v_hllc0;
        }

    for (int iy = 0; iy < OutputSOA::NY; ++iy)
        for (int ix = 0; ix < OutputSOA::NX; ++ix)
        {
            const Real a2p = A2p(iy,ix);
            const Real a2m = A2m(iy+1,ix);
            const Real a1p = static_cast<Real>(1.0) - a2p;
            const Real a1m = static_cast<Real>(1.0) - a2m;

            // A2
            //sumA2.ref(ix, iy) += A2p(iy, ix) + A2m(iy+1, ix);
            sumA2.ref(ix, iy) += a2p + a2m;

            // K
#ifndef LAMBDA
            // TODO: (fabianw@mavt.ethz.ch; Mon 02 May 2016 09:26:04 AM CEST)
            // testing
            const Real num1m = max(static_cast<Real>(0.0), _g1*(pm(iy+1,ix) + _pc1));
            const Real num2m = max(static_cast<Real>(0.0), _g2*(pm(iy+1,ix) + _pc2));
            const Real num1p = max(static_cast<Real>(0.0), _g1*(pp(iy,ix) + _pc1));
            const Real num2p = max(static_cast<Real>(0.0), _g2*(pp(iy,ix) + _pc2));
            // sumK.ref(ix,iy) += (0 < a1m && a1m < 1) ? (num1m - num2m)/(num1m/a1m + num2m/a2m) : 0;
            // sumK.ref(ix,iy) += (0 < a1p && a1p < 1) ? (num1p - num2p)/(num1p/a1p + num2p/a2p) : 0;
            sumK.ref(ix,iy) += (static_cast<Real>(ALPHAEPS) < a1m && a1m < static_cast<Real>(1-ALPHAEPS)) ? (num1m - num2m)/(num1m/a1m + num2m/a2m) : 0;
            sumK.ref(ix,iy) += (static_cast<Real>(ALPHAEPS) < a1p && a1p < static_cast<Real>(1-ALPHAEPS)) ? (num1p - num2p)/(num1p/a1p + num2p/a2p) : 0;

            // sumK.ref(ix,iy) += a1p*a2p*(pp(iy,ix)*(_g1-_g2) + _g1*_pc1 - _g2*_pc2)/(pp(iy,ix)*(a1p*_g2 + a2p*_g1) + a1p*_g2*_pc2 + a2p*_g1*_pc1)
            // + a1m*a2m*(pm(iy+1,ix)*(_g1-_g2) + _g1*_pc1 - _g2*_pc2)/(pm(iy+1,ix)*(a1m*_g2 + a2m*_g1) + a1m*_g2*_pc2 + a2m*_g1*_pc1);
#else
            const Real gpc1m = _g1*(pm(iy+1,ix) + _pc1);
            const Real gpc2m = _g2*(pm(iy+1,ix) + _pc2);
            const Real gpc1p = _g1*(pp(iy,ix) + _pc1);
            const Real gpc2p = _g2*(pp(iy,ix) + _pc2);

            const Real sela2m = (a2m > 0) ? 1.0 : -1.0;
            const Real sela2p = (a2p > 0) ? 1.0 : -1.0;
            const Real selgpc1m = (gpc1m > 0) ? 1.0 : -1.0;
            const Real selgpc2m = (gpc2m > 0) ? 1.0 : -1.0;
            const Real selgpc1p = (gpc1p > 0) ? 1.0 : -1.0;
            const Real selgpc2p = (gpc2p > 0) ? 1.0 : -1.0;

            const Real a2m_pos   = (sela2m < 0)   ? 1.0 : a2m;
            const Real gpc1m_pos = (selgpc1m < 0) ? 1.0 : gpc1m;
            const Real gpc2m_pos = (selgpc2m < 0) ? 0.0 : gpc2m;
            const Real a2p_pos   = (sela2p < 0)   ? 1.0 : a2p;
            const Real gpc1p_pos = (selgpc1p < 0) ? 1.0 : gpc1p;
            const Real gpc2p_pos = (selgpc2p < 0) ? 0.0 : gpc2p;

            const Real a1m_pos = 1.0 - a2m_pos;
            const Real a1p_pos = 1.0 - a2p_pos;
            const Real lambdam = 0.25*(1.0+sela2m)*(1.0+selgpc1m);
            const Real lambdap = 0.25*(1.0+sela2p)*(1.0+selgpc1p);

            sumK.ref(ix,iy) += lambdam*a1m_pos*a2m_pos*(gpc1m_pos - gpc2m_pos) / (a1m_pos*gpc2m_pos + a2m_pos*gpc1m_pos);
            sumK.ref(ix,iy) += lambdap*a1p_pos*a2p_pos*(gpc1p_pos - gpc2p_pos) / (a1p_pos*gpc2p_pos + a2p_pos*gpc1p_pos);
#endif /* 0 */
        }
}

void Convection_CPP_HLLC_5eq::_zextraterm_hllc( const TempSOA& wm0, const TempSOA& wp0,
        const TempSOA& wm1, const TempSOA& wp1,
        const TempSOA& A2m,  const TempSOA& A2p,
        const TempSOA& pm, const TempSOA& pp,
        const TempSOA& sm0, const TempSOA& sp0, const TempSOA& ss0,
        const TempSOA& sm1, const TempSOA& sp1, const TempSOA& ss1)
{
    for (int iy = 0; iy < OutputSOA::NY; ++iy)
        for (int ix = 0; ix < OutputSOA::NX; ++ix)
        {
            int sign_star;
            Real s_minus, s_pluss, chi_starm, chi_starp;

            // right cell interface (w_{i+1/2})
            sign_star = ss1(ix, iy) == 0 ? 0 : (ss1(ix, iy) < 0 ? -1 : 1);
            s_minus   = min(static_cast<Real>(0.0), sm1(ix, iy));
            s_pluss   = max(static_cast<Real>(0.0), sp1(ix, iy));
            chi_starm = (sm1(ix, iy) - wm1(ix, iy))/(sm1(ix, iy) - ss1(ix, iy)) - static_cast<Real>(1.0);
            chi_starp = (sp1(ix, iy) - wp1(ix, iy))/(sp1(ix, iy) - ss1(ix, iy)) - static_cast<Real>(1.0);
            const Real w_hllc1 = static_cast<Real>(0.5*(1 + sign_star))*(wm1(ix, iy) + s_minus*chi_starm) + static_cast<Real>(0.5*(1 - sign_star))*(wp1(ix, iy) + s_pluss*chi_starp);

            // left cell interface (w_{i-1/2})
            sign_star = ss0(ix, iy) == 0 ? 0 : (ss0(ix, iy) < 0 ? -1 : 1);
            s_minus   = min(static_cast<Real>(0.0), sm0(ix, iy));
            s_pluss   = max(static_cast<Real>(0.0), sp0(ix, iy));
            chi_starm = (sm0(ix, iy) - wm0(ix, iy))/(sm0(ix, iy) - ss0(ix, iy)) - static_cast<Real>(1.0);
            chi_starp = (sp0(ix, iy) - wp0(ix, iy))/(sp0(ix, iy) - ss0(ix, iy)) - static_cast<Real>(1.0);
            const Real w_hllc0 = static_cast<Real>(0.5*(1 + sign_star))*(wm0(ix, iy) + s_minus*chi_starm) + static_cast<Real>(0.5*(1 - sign_star))*(wp0(ix, iy) + s_pluss*chi_starp);

            divu.ref(ix, iy) += w_hllc1 - w_hllc0;
        }

    for (int iy = 0; iy < OutputSOA::NY; ++iy)
        for (int ix = 0; ix < OutputSOA::NX; ++ix)
        {
            const Real a2p = A2p(ix,iy);
            const Real a2m = A2m(ix,iy);
            const Real a1p = static_cast<Real>(1.0) - a2p;
            const Real a1m = static_cast<Real>(1.0) - a2m;

            // A2
            //sumA2.ref(ix, iy) += A2p(ix, iy) + A2m(ix, iy);
            sumA2.ref(ix, iy) += a2p + a2m;

            // K
#ifndef LAMBDA
            // TODO: (fabianw@mavt.ethz.ch; Mon 02 May 2016 09:26:04 AM CEST)
            // testing
            const Real num1m = max(static_cast<Real>(0.0), _g1*(pm(ix,iy) + _pc1));
            const Real num2m = max(static_cast<Real>(0.0), _g2*(pm(ix,iy) + _pc2));
            const Real num1p = max(static_cast<Real>(0.0), _g1*(pp(ix,iy) + _pc1));
            const Real num2p = max(static_cast<Real>(0.0), _g2*(pp(ix,iy) + _pc2));
            // sumK.ref(ix,iy) += (0 < a1m && a1m < 1) ? (num1m - num2m)/(num1m/a1m + num2m/a2m) : 0;
            // sumK.ref(ix,iy) += (0 < a1p && a1p < 1) ? (num1p - num2p)/(num1p/a1p + num2p/a2p) : 0;
            sumK.ref(ix,iy) += (static_cast<Real>(ALPHAEPS) < a1m && a1m < static_cast<Real>(1-ALPHAEPS)) ? (num1m - num2m)/(num1m/a1m + num2m/a2m) : 0;
            sumK.ref(ix,iy) += (static_cast<Real>(ALPHAEPS) < a1p && a1p < static_cast<Real>(1-ALPHAEPS)) ? (num1p - num2p)/(num1p/a1p + num2p/a2p) : 0;

            // sumK.ref(ix,iy) += a1p*a2p*(pp(ix,iy)*(_g1-_g2) + _g1*_pc1 - _g2*_pc2)/(pp(ix,iy)*(a1p*_g2 + a2p*_g1) + a1p*_g2*_pc2 + a2p*_g1*_pc1)
            // + a1m*a2m*(pm(ix,iy)*(_g1-_g2) + _g1*_pc1 - _g2*_pc2)/(pm(ix,iy)*(a1m*_g2 + a2m*_g1) + a1m*_g2*_pc2 + a2m*_g1*_pc1);
#else
            const Real gpc1m = _g1*(pm(ix,iy) + _pc1);
            const Real gpc2m = _g2*(pm(ix,iy) + _pc2);
            const Real gpc1p = _g1*(pp(ix,iy) + _pc1);
            const Real gpc2p = _g2*(pp(ix,iy) + _pc2);

            const Real sela2m = (a2m > 0) ? 1.0 : -1.0;
            const Real sela2p = (a2p > 0) ? 1.0 : -1.0;
            const Real selgpc1m = (gpc1m > 0) ? 1.0 : -1.0;
            const Real selgpc2m = (gpc2m > 0) ? 1.0 : -1.0;
            const Real selgpc1p = (gpc1p > 0) ? 1.0 : -1.0;
            const Real selgpc2p = (gpc2p > 0) ? 1.0 : -1.0;

            const Real a2m_pos   = (sela2m < 0)   ? 1.0 : a2m;
            const Real gpc1m_pos = (selgpc1m < 0) ? 1.0 : gpc1m;
            const Real gpc2m_pos = (selgpc2m < 0) ? 0.0 : gpc2m;
            const Real a2p_pos   = (sela2p < 0)   ? 1.0 : a2p;
            const Real gpc1p_pos = (selgpc1p < 0) ? 1.0 : gpc1p;
            const Real gpc2p_pos = (selgpc2p < 0) ? 0.0 : gpc2p;

            const Real a1m_pos = 1.0 - a2m_pos;
            const Real a1p_pos = 1.0 - a2p_pos;
            const Real lambdam = 0.25*(1.0+sela2m)*(1.0+selgpc1m);
            const Real lambdap = 0.25*(1.0+sela2p)*(1.0+selgpc1p);

            sumK.ref(ix,iy) += lambdam*a1m_pos*a2m_pos*(gpc1m_pos - gpc2m_pos) / (a1m_pos*gpc2m_pos + a2m_pos*gpc1m_pos);
            sumK.ref(ix,iy) += lambdap*a1p_pos*a2p_pos*(gpc1p_pos - gpc2p_pos) / (a1p_pos*gpc2p_pos + a2p_pos*gpc1p_pos);
#endif /* 0 */
        }
}


// ============================================================================
// Characteristic velocity based on Einfeldts work (1988)
void Convection_CPP_HLLC_5eq::_char_vel_einfeldt( const int SizeX,
        const TempSOA& r1m, const TempSOA& r1p,
        const TempSOA& r2m, const TempSOA& r2p,
        const TempSOA& vm, const TempSOA& vp,
        const TempSOA& pm, const TempSOA& pp,
        const TempSOA& A2m, const TempSOA& A2p,
        TempSOA& outm, TempSOA& outp) // ? FLOP
{
    /* *
     * Compute upper and lower bounds of signal velocities for the Riemann
     * problem according to Einfeldt:
     *
     * 1.) Compute Rr needed for Roe averages
     * 2.) Compute speed of sound in left and right state
     * 3.) Compute speed of sound according to Einfeldt and Rr
     * 4.) Compute upper and lower signal velocities
     * */
    for (int iy = 0; iy < TempSOA::NY; ++iy)
        for (int ix = 0; ix < SizeX; ++ix)
        {
            // 1.)
            const Real rm = r1m(ix,iy) + r2m(ix,iy);
            const Real rp = r1p(ix,iy) + r2p(ix,iy);
            const Real Rr   = mysqrt(rp/rm);
            const Real Rinv = 1.0 / (1.0 + Rr);

            // 2.)
            const Real alpha2m = A2m(ix,iy);
            const Real alpha2p = A2p(ix,iy);

            const Real alpha1m = static_cast<Real>(1.0)-alpha2m;
            const Real alpha1p = static_cast<Real>(1.0)-alpha2p;

#ifdef _NOK_
            const Real GmixInvm = static_cast<Real>(1.0) / (alpha1m*_g1m1Inv + alpha2m*_g2m1Inv);
            const Real GmixInvp = static_cast<Real>(1.0) / (alpha1p*_g1m1Inv + alpha2p*_g2m1Inv);
            const Real Pmixm = alpha1m*_g1m1Inv*_g1*_pc1 + alpha2m*_g2m1Inv*_g2*_pc2;
            const Real Pmixp = alpha1p*_g1m1Inv*_g1*_pc1 + alpha2p*_g2m1Inv*_g2*_pc2;

            const Real cm2 = ((GmixInvm + static_cast<Real>(1.0))*pm(ix,iy) + GmixInvm*Pmixm)/rm;
            const Real cp2 = ((GmixInvp + static_cast<Real>(1.0))*pp(ix,iy) + GmixInvp*Pmixp)/rp;
#else

            // TODO: (fabianw@mavt.ethz.ch; Fri 29 Apr 2016 06:25:35 PM CEST)
            // testing different numeric computation of mixture speed of sound
            const Real denom1_m = _g1*(pm(ix,iy) + _pc1);
            const Real denom2_m = _g2*(pm(ix,iy) + _pc2);
            const Real denom1_p = _g1*(pp(ix,iy) + _pc1);
            const Real denom2_p = _g2*(pp(ix,iy) + _pc2);
            Real rc2m = (denom1_m > 0) ? alpha1m/denom1_m : 0;
            rc2m += (denom2_m > 0) ? alpha2m/denom2_m : 0;
            Real rc2p = (denom1_p > 0) ? alpha1p/denom1_p : 0;
            rc2p += (denom2_p > 0) ? alpha2p/denom2_p : 0;
            const Real cm2 = 1.0/(rm*rc2m);
            const Real cp2 = 1.0/(rp*rc2p);

            // const Real rc2m = _g1*_g2*(pm(ix,iy) + _pc1)*(pm(ix,iy) + _pc2)/(alpha1m*_g2*(pm(ix,iy) + _pc2) + alpha2m*_g1*(pm(ix,iy) + _pc1));
            // const Real rc2p = _g1*_g2*(pp(ix,iy) + _pc1)*(pp(ix,iy) + _pc2)/(alpha1p*_g2*(pp(ix,iy) + _pc2) + alpha2p*_g1*(pp(ix,iy) + _pc1));
            // const Real cm2 = rc2m/rm;
            // const Real cp2 = rc2p/rp;
#endif /* _NOK_ */

            const Real cm = mysqrt(cm2);
            const Real cp = mysqrt(cp2);

#ifndef NDEBUG
#ifdef _NOK_
            if(std::isnan(cm))
                //printf("cm hllc fails now. Ingredients are %e %e %e %e %e\n", pm(ix,iy), GmixInvm, Pmixm, rm, alpha2m);
                cout << "cm hllc fails now. Ingredients are " << setprecision(12) << pm(ix,iy) << " " << GmixInvm << " " << Pmixm << " " << rm << " " << alpha2m << endl;

            if(std::isnan(cp))
                //printf("cp hllc fails now. Ingredients are %e %e %e %e %e\n", pp(ix,iy), GmixInvp, Pmixp, rp, alpha2p);
                cout << "cp hllc fails now. Ingredients are " << setprecision(12) << pp(ix,iy) << " " << GmixInvp << " " << Pmixp << " " << rp << " " << alpha2p << endl;
#else
            if(std::isnan(cm))
                printf("cm hllc fails now. Ingredients are %e %e %e %e\n", pm(ix,iy), rc2m, rm, alpha2m);
            if(std::isnan(cp))
                printf("cp hllc fails now. Ingredients are %e %e %e %e\n", pp(ix,iy), rc2p, rp, alpha2p);
#endif /* _NOK_ */
#endif

            assert(!std::isnan(cm));
            assert(!std::isnan(cp));

            // 3.)
            const Real um    = vm(ix, iy);
            const Real up    = vp(ix, iy);
            const Real eta_2 = 0.5*Rr*Rinv*Rinv;
            const Real d2    = (cm2 + Rr*cp2)*Rinv + eta_2*(up - um)*(up - um);
            const Real d     = mysqrt(d2);
            const Real u     = (um + Rr*up)*Rinv;

            assert(!std::isnan(d));
            assert(!std::isnan(u));

            // 4.)
            outm.ref(ix, iy) = min(u - d, um - cm);
            outp.ref(ix, iy) = max(u + d, up + cp);
        }
}


/* *
 * Compute characteristic velocity, s^star, of the intermediate wave.  The
 * computation is based on the condition of uniform constant pressure in
 * the star region.  See P. Batten et. al., "On the choice of wavespeeds
 * for the HLLC Riemann solver", SIAM J. Sci. Comput. 18 (1997) 1553--1570
 * It is assumed s^minus and s^plus are known.
 * */
void Convection_CPP_HLLC_5eq::_char_vel_star(const int SizeX,
        const TempSOA& r1m, const TempSOA& r1p,
        const TempSOA& r2m, const TempSOA& r2p,
        const TempSOA& vm, const TempSOA& vp,
        const TempSOA& pm, const TempSOA& pp,
        const TempSOA& sm, const TempSOA& sp,
        TempSOA& out_star) // 11 FLOP
{
    for (int iy = 0; iy < TempSOA::NY; ++iy)
        for (int ix = 0; ix < SizeX; ++ix)
        {
            const Real facm = (r1m(ix, iy)+r2m(ix, iy)) * (sm(ix, iy) - vm(ix, iy));
            const Real facp = (r1p(ix, iy)+r2p(ix, iy)) * (sp(ix, iy) - vp(ix, iy));
            out_star.ref(ix, iy) = (pp(ix, iy) - pm(ix, iy) + vm(ix, iy)*facm - vp(ix, iy)*facp)/(facm - facp);

            assert(!std::isnan(out_star.ref(ix, iy)));
        }
}

// ============================================================================


// HLLC FLUXES
void Convection_CPP_HLLC_5eq::_hllc_rho( const int SizeX,
        const TempSOA& rm, const TempSOA& rp,
        const TempSOA& vm, const TempSOA& vp,
        const TempSOA& sm, const TempSOA& sp, const TempSOA& ss,
        TempSOA& flux_out) // 23 FLOP
{
    for (int iy = 0; iy < TempSOA::NY; ++iy)
        for (int ix = 0; ix < SizeX; ++ix)
        {
            /* *
             * The flux computation is split into 4 parts:
             * 1.) Compute signum of s^*, compute s^- and s^+
             * 2.) Compute chi^* and delta of q^* and q
             * 3.) Compute trivial flux
             * 4.) Compute HLLC flux
             * */

            // 1.)
            const int sign_star = ss(ix, iy) == 0 ? 0 : (ss(ix, iy) < 0 ? -1 : 1);
            const Real s_minus  = min(static_cast<Real>(0.0), sm(ix, iy));
            const Real s_pluss  = max(static_cast<Real>(0.0), sp(ix, iy));

            // 2.)
            const Real chi_starm = (sm(ix, iy) - vm(ix, iy)) / (sm(ix, iy) - ss(ix, iy));
            const Real chi_starp = (sp(ix, iy) - vp(ix, iy)) / (sp(ix, iy) - ss(ix, iy));
            const Real qm        = rm(ix, iy);
            const Real qp        = rp(ix, iy);
            const Real q_deltam  = qm*chi_starm - qm;
            const Real q_deltap  = qp*chi_starp - qp;

            // 3.)
            const Real fm = qm*vm(ix, iy);
            const Real fp = qp*vp(ix, iy);

            // 4.)
            flux_out.ref(ix, iy) = static_cast<Real>(0.5*(1 + sign_star))*(fm + s_minus*q_deltam) + static_cast<Real>(0.5*(1 - sign_star))*(fp + s_pluss*q_deltap);
            assert(!std::isnan(flux_out.ref(ix, iy)));
        }
}

void Convection_CPP_HLLC_5eq::_hllc_vel(const int SizeX,
        const TempSOA& r1m,  const TempSOA& r1p,
        const TempSOA& r2m,  const TempSOA& r2p,
        const TempSOA& vm,  const TempSOA& vp,
        const TempSOA& vdm, const TempSOA& vdp,
        const TempSOA& sm,  const TempSOA& sp,  const TempSOA& ss,
        TempSOA& flux_out) // 25 FLOP
{
    for (int iy = 0; iy < TempSOA::NY; ++iy)
        for (int ix = 0; ix < SizeX; ++ix)
        {
            /* *
             * The flux computation is split into 4 parts:
             * 1.) Compute signum of s^*, compute s^- and s^+
             * 2.) Compute chi^* and delta of q^* and q
             * 3.) Compute trivial flux
             * 4.) Compute HLLC flux
             * */

            // 1.)
            const int sign_star = ss(ix, iy) == 0 ? 0 : (ss(ix, iy) < 0 ? -1 : 1);
            const Real s_minus  = min(static_cast<Real>(0.0), sm(ix, iy));
            const Real s_pluss  = max(static_cast<Real>(0.0), sp(ix, iy));

            // 2.)
            const Real rm = r1m(ix,iy) + r2m(ix,iy);
            const Real rp = r1p(ix,iy) + r2p(ix,iy);
            const Real chi_starm = (sm(ix, iy) - vdm(ix, iy)) / (sm(ix, iy) - ss(ix, iy));
            const Real chi_starp = (sp(ix, iy) - vdp(ix, iy)) / (sp(ix, iy) - ss(ix, iy));
            const Real qm        = rm*vm(ix, iy);
            const Real qp        = rp*vp(ix, iy);
            const Real q_deltam  = qm*chi_starm - qm;
            const Real q_deltap  = qp*chi_starp - qp;

            // 3.)
            const Real fm = qm*vdm(ix, iy);
            const Real fp = qp*vdp(ix, iy);

            // 4.)
            flux_out.ref(ix, iy) = static_cast<Real>(0.5*(1 + sign_star))*(fm + s_minus*q_deltam) + static_cast<Real>(0.5*(1 - sign_star))*(fp + s_pluss*q_deltap);
            assert(!std::isnan(flux_out.ref(ix, iy)));
            assert(!std::isnan(ss(ix, iy)));
            assert(!std::isnan(sm(ix, iy)));
            assert(!std::isnan(sp(ix, iy)));
        }
}

void Convection_CPP_HLLC_5eq::_hllc_pvel(const int SizeX,
        const TempSOA& r1m,  const TempSOA& r1p,
        const TempSOA& r2m,  const TempSOA& r2p,
        const TempSOA& vm, const TempSOA& vp,
        const TempSOA& pm, const TempSOA& pp,
        const TempSOA& sm, const TempSOA& sp, const TempSOA& ss,
        TempSOA& flux_out) // 29 FLOP
{
    for (int iy = 0; iy < TempSOA::NY; ++iy)
        for (int ix = 0; ix < SizeX; ++ix)
        {
            /* *
             * The flux computation is split into 4 parts:
             * 1.) Compute signum of s^*, compute s^- and s^+
             * 2.) Compute chi^* and delta of q^* and q
             * 3.) Compute trivial flux
             * 4.) Compute HLLC flux
             * */

            // 1.)
            const int sign_star = ss(ix, iy) == 0 ? 0 : (ss(ix, iy) < 0 ? -1 : 1);
            const Real s_minus  = min(static_cast<Real>(0.0), sm(ix, iy));
            const Real s_pluss  = max(static_cast<Real>(0.0), sp(ix, iy));

            // 2.)
            const Real rm = r1m(ix,iy) + r2m(ix,iy);
            const Real rp = r1p(ix,iy) + r2p(ix,iy);
            const Real chi_starm = (sm(ix, iy) - vm(ix, iy)) / (sm(ix, iy) - ss(ix, iy));
            const Real chi_starp = (sp(ix, iy) - vp(ix, iy)) / (sp(ix, iy) - ss(ix, iy));
            const Real qm        = rm*vm(ix, iy);
            const Real qp        = rp*vp(ix, iy);
            const Real q_deltam  = rm*ss(ix, iy)*chi_starm - qm;
            const Real q_deltap  = rp*ss(ix, iy)*chi_starp - qp;

            // 3.)
            const Real fm = qm*vm(ix, iy) + pm(ix, iy);
            const Real fp = qp*vp(ix, iy) + pp(ix, iy);

            // 4.)
            flux_out.ref(ix, iy) = static_cast<Real>(0.5*(1 + sign_star))*(fm + s_minus*q_deltam) + static_cast<Real>(0.5*(1 - sign_star))*(fp + s_pluss*q_deltap);
            assert(!std::isnan(flux_out.ref(ix, iy)));
            assert(rm > 0);
            assert(rp > 0);
        }
}

void Convection_CPP_HLLC_5eq::_hllc_e(const int SizeX,
        const TempSOA& r1m, const TempSOA& r1p,
        const TempSOA& r2m, const TempSOA& r2p,
        const TempSOA& vdm, const TempSOA& vdp,
        const TempSOA& v1m, const TempSOA& v1p,
        const TempSOA& v2m, const TempSOA& v2p,
        const TempSOA& pm,  const TempSOA& pp,
        const TempSOA& A2m, const TempSOA& A2p,
        const TempSOA& sm,  const TempSOA& sp,  const TempSOA& ss,
        TempSOA& flux_out) // 59 FLOP
{
    for (int iy = 0; iy < TempSOA::NY; ++iy)
        for (int ix = 0; ix < SizeX; ++ix)
        {
            /* *
             * The flux computation is split into 4 parts:
             * 1.) Compute signum of s^*, compute s^- and s^+
             * 2.) Compute helpers
             * 3.) Compute chi^* and delta of q^* and q
             * 4.) Compute trivial flux
             * 5.) Compute HLLC flux
             * */

            // 1.)
            const int sign_star = ss(ix, iy) == 0 ? 0 : (ss(ix, iy) < 0 ? -1 : 1);
            const Real alpha2m  = A2m(ix,iy);
            const Real alpha2p  = A2p(ix,iy);
            const Real alpha1m  = static_cast<Real>(1.0) - alpha2m;
            const Real alpha1p  = static_cast<Real>(1.0) - alpha2p;
            const Real s_minus  = min(static_cast<Real>(0.0), sm(ix, iy));
            const Real s_pluss  = max(static_cast<Real>(0.0), sp(ix, iy));

            // 2.)
            const Real rm = r1m(ix,iy) + r2m(ix,iy);
            const Real rp = r1p(ix,iy) + r2p(ix,iy);
            const Real gmix_m1Invm = alpha1m * _g1m1Inv + alpha2m * _g2m1Inv;
            const Real gmix_m1Invp = alpha1p * _g1m1Inv + alpha2p * _g2m1Inv;
            const Real pcm = max(static_cast<Real>(0.0),alpha1m *  _g1m1Inv*_g1*_pc1 + alpha2m *  _g2m1Inv*_g2*_pc2);
            const Real pcp = max(static_cast<Real>(0.0),alpha1p *  _g1m1Inv*_g1*_pc1 + alpha2p *  _g2m1Inv*_g2*_pc2);

            // 3.)
            const Real chi_starm = (sm(ix, iy) - vdm(ix, iy)) / (sm(ix, iy) - ss(ix, iy));
            const Real chi_starp = (sp(ix, iy) - vdp(ix, iy)) / (sp(ix, iy) - ss(ix, iy));
            const Real qm        = gmix_m1Invm*pm(ix, iy) + pcm + static_cast<Real>(0.5)*rm*(vdm(ix, iy)*vdm(ix, iy) + v1m(ix, iy)*v1m(ix, iy) + v2m(ix, iy)*v2m(ix, iy));
            const Real qp        = gmix_m1Invp*pp(ix, iy) + pcp + static_cast<Real>(0.5)*rp*(vdp(ix, iy)*vdp(ix, iy) + v1p(ix, iy)*v1p(ix, iy) + v2p(ix, iy)*v2p(ix, iy));
            const Real q_deltam  = chi_starm*(qm + (ss(ix, iy) - vdm(ix, iy))*(rm*ss(ix, iy) + pm(ix, iy)/(sm(ix, iy) - vdm(ix, iy)))) - qm;
            const Real q_deltap  = chi_starp*(qp + (ss(ix, iy) - vdp(ix, iy))*(rp*ss(ix, iy) + pp(ix, iy)/(sp(ix, iy) - vdp(ix, iy)))) - qp;

            // 4.)
            const Real fm = vdm(ix, iy)*(qm + pm(ix, iy));
            const Real fp = vdp(ix, iy)*(qp + pp(ix, iy));

            // 5.)
            flux_out.ref(ix, iy) = static_cast<Real>(0.5*(1 + sign_star))*(fm + s_minus*q_deltam) + static_cast<Real>(0.5*(1 - sign_star))*(fp + s_pluss*q_deltap);
            assert(!std::isnan(flux_out.ref(ix, iy)));
        }
}


// VIRTUAL METHODS
void Convection_CPP_HLLC_5eq::_xflux(const int relid)
{
#ifndef _A2_A2R2_R_
    _xweno_minus_clipped(rho1.ring(relid), rho1.weno.ref(0));
    _xweno_pluss_clipped(rho1.ring(relid), rho1.weno.ref(1));
#endif
    _xweno_minus_clipped(rho2.ring(relid), rho2.weno.ref(0));
    _xweno_pluss_clipped(rho2.ring(relid), rho2.weno.ref(1));
#ifdef _A2_A2R2_R_
    _xweno_minus_clipped(rho.ring(relid), rho.weno.ref(0));
    _xweno_pluss_clipped(rho.ring(relid), rho.weno.ref(1));
#endif
    _xweno_minus_clipped(u.ring(relid), u.weno.ref(0));
    _xweno_pluss_clipped(u.ring(relid), u.weno.ref(1));
    _xweno_minus_clipped(v.ring(relid), v.weno.ref(0));
    _xweno_pluss_clipped(v.ring(relid), v.weno.ref(1));
    _xweno_minus_clipped(w.ring(relid), w.weno.ref(0));
    _xweno_pluss_clipped(w.ring(relid), w.weno.ref(1));
    _xweno_minus_clipped(p.ring(relid), p.weno.ref(0));
    _xweno_pluss_clipped(p.ring(relid), p.weno.ref(1));
    _xweno_minus_clipped(A2.ring(relid), A2.weno.ref(0));
    _xweno_pluss_clipped(A2.ring(relid), A2.weno.ref(1));

#ifdef _A2_A2R2_R_
    _convert_reconstruction();
#endif

#ifdef _RECONPCLIP_
   _pressure_clipping_hllc(TempSOA::NX, A2.weno(0), A2.weno(1),p.weno.ref(0), p.weno.ref(1));
#endif

    // characteristic velocities
    _char_vel_einfeldt(TempSOA::NX, rho1.weno(0), rho1.weno(1), rho2.weno(0), rho2.weno(1), u.weno(0), u.weno(1), p.weno(0), p.weno(1), A2.weno(0), A2.weno(1), charvel.ref(0), charvel.ref(1));

    // intermediate wave
    _char_vel_star(TempSOA::NX, rho1.weno(0), rho1.weno(1), rho2.weno(0), rho2.weno(1), u.weno(0), u.weno(1), p.weno(0), p.weno(1), charvel(0), charvel(1), charvel_star.ref(0));

    // RHS of advection equations
    /* _xextraterm_hllc(u.weno(0), u.weno(1), A2.weno(0), A2.weno(1), charvel(0), charvel(1), charvel_star(0)); */
    _xextraterm_hllc(u.weno(0), u.weno(1), A2.weno(0), A2.weno(1), p.weno(0), p.weno(1), charvel(0), charvel(1), charvel_star(0));

    // conservative variables
    _hllc_rho( TempSOA::NX, rho1.weno(0), rho1.weno(1), u.weno(0), u.weno(1), charvel(0), charvel(1), charvel_star(0), rho1.flux.ref());
    _hllc_rho( TempSOA::NX, rho2.weno(0), rho2.weno(1), u.weno(0), u.weno(1), charvel(0), charvel(1), charvel_star(0), rho2.flux.ref());
    _hllc_pvel(TempSOA::NX, rho1.weno(0), rho1.weno(1), rho2.weno(0), rho2.weno(1), u.weno(0), u.weno(1), p.weno(0), p.weno(1), charvel(0), charvel(1), charvel_star(0), u.flux.ref());
    _hllc_vel( TempSOA::NX, rho1.weno(0), rho1.weno(1), rho2.weno(0), rho2.weno(1), v.weno(0), v.weno(1), u.weno(0), u.weno(1), charvel(0), charvel(1), charvel_star(0), v.flux.ref());
    _hllc_vel( TempSOA::NX, rho1.weno(0), rho1.weno(1), rho2.weno(0), rho2.weno(1), w.weno(0), w.weno(1), u.weno(0), u.weno(1), charvel(0), charvel(1), charvel_star(0), w.flux.ref());
    _hllc_e(   TempSOA::NX, rho1.weno(0), rho1.weno(1), rho2.weno(0), rho2.weno(1), u.weno(0), u.weno(1), v.weno(0), v.weno(1), w.weno(0), w.weno(1), p.weno(0), p.weno(1), A2.weno(0), A2.weno(1), charvel(0), charvel(1), charvel_star(0), p.flux.ref());

    // advected quantities
    _hllc_rho(TempSOA::NX, A2.weno(0), A2.weno(1), u.weno(0), u.weno(1), charvel(0), charvel(1), charvel_star(0), A2.flux.ref());
}

void Convection_CPP_HLLC_5eq::_yflux(const int relid)
{
#ifndef _A2_A2R2_R_
    _yweno_minus_clipped(rho1.ring(relid), rho1.weno.ref(0));
    _yweno_pluss_clipped(rho1.ring(relid), rho1.weno.ref(1));
#endif
    _yweno_minus_clipped(rho2.ring(relid), rho2.weno.ref(0));
    _yweno_pluss_clipped(rho2.ring(relid), rho2.weno.ref(1));
#ifdef _A2_A2R2_R_
    _yweno_minus_clipped(rho.ring(relid), rho.weno.ref(0));
    _yweno_pluss_clipped(rho.ring(relid), rho.weno.ref(1));
#endif
    _yweno_minus_clipped(u.ring(relid), u.weno.ref(0));
    _yweno_pluss_clipped(u.ring(relid), u.weno.ref(1));
    _yweno_minus_clipped(v.ring(relid), v.weno.ref(0));
    _yweno_pluss_clipped(v.ring(relid), v.weno.ref(1));
    _yweno_minus_clipped(w.ring(relid), w.weno.ref(0));
    _yweno_pluss_clipped(w.ring(relid), w.weno.ref(1));
    _yweno_minus_clipped(p.ring(relid), p.weno.ref(0));
    _yweno_pluss_clipped(p.ring(relid), p.weno.ref(1));
    _yweno_minus_clipped(A2.ring(relid), A2.weno.ref(0));
    _yweno_pluss_clipped(A2.ring(relid), A2.weno.ref(1));

#ifdef _A2_A2R2_R_
    _convert_reconstruction();
#endif

#ifdef _RECONPCLIP_
   _pressure_clipping_hllc(TempSOA::NX, A2.weno(0), A2.weno(1),p.weno.ref(0), p.weno.ref(1));
#endif

    // characteristic velocities
    _char_vel_einfeldt(TempSOA::NX, rho1.weno(0), rho1.weno(1), rho2.weno(0), rho2.weno(1), v.weno(0), v.weno(1), p.weno(0), p.weno(1), A2.weno(0), A2.weno(1), charvel.ref(0), charvel.ref(1));

    // intermediate wave
    _char_vel_star(TempSOA::NX, rho1.weno(0), rho1.weno(1), rho2.weno(0), rho2.weno(1), v.weno(0), v.weno(1), p.weno(0), p.weno(1), charvel(0), charvel(1), charvel_star.ref(0));

    // RHS of advection equations
    /* _yextraterm_hllc(v.weno(0), v.weno(1), A2.weno(0), A2.weno(1), charvel(0), charvel(1), charvel_star(0)); */
    _yextraterm_hllc(v.weno(0), v.weno(1), A2.weno(0), A2.weno(1), p.weno(0), p.weno(1), charvel(0), charvel(1), charvel_star(0));

    // conservative variables
    _hllc_rho( TempSOA::NX, rho1.weno(0), rho1.weno(1), v.weno(0), v.weno(1), charvel(0), charvel(1), charvel_star(0), rho1.flux.ref());
    _hllc_rho( TempSOA::NX, rho2.weno(0), rho2.weno(1), v.weno(0), v.weno(1), charvel(0), charvel(1), charvel_star(0), rho2.flux.ref());
    _hllc_vel( TempSOA::NX, rho1.weno(0), rho1.weno(1), rho2.weno(0), rho2.weno(1), u.weno(0), u.weno(1), v.weno(0), v.weno(1), charvel(0), charvel(1), charvel_star(0), u.flux.ref());
    _hllc_pvel(TempSOA::NX, rho1.weno(0), rho1.weno(1), rho2.weno(0), rho2.weno(1), v.weno(0), v.weno(1), p.weno(0), p.weno(1), charvel(0), charvel(1), charvel_star(0), v.flux.ref());
    _hllc_vel( TempSOA::NX, rho1.weno(0), rho1.weno(1), rho2.weno(0), rho2.weno(1), w.weno(0), w.weno(1), v.weno(0), v.weno(1), charvel(0), charvel(1), charvel_star(0), w.flux.ref());
    _hllc_e(   TempSOA::NX, rho1.weno(0), rho1.weno(1), rho2.weno(0), rho2.weno(1), v.weno(0), v.weno(1), u.weno(0), u.weno(1), w.weno(0), w.weno(1), p.weno(0), p.weno(1), A2.weno(0), A2.weno(1), charvel(0), charvel(1), charvel_star(0), p.flux.ref());

    // advected quantities
    _hllc_rho(TempSOA::NX, A2.weno(0), A2.weno(1), v.weno(0), v.weno(1), charvel(0), charvel(1), charvel_star(0), A2.flux.ref());
}

void Convection_CPP_HLLC_5eq::_zflux(const int relid, const bool evalXtra)
{
#ifndef _A2_A2R2_R_
    _zweno_minus_clipped(relid, rho1.ring, rho1.weno.ref(0));
    _zweno_pluss_clipped(relid, rho1.ring, rho1.weno.ref(1));
#endif
    _zweno_minus_clipped(relid, rho2.ring, rho2.weno.ref(0));
    _zweno_pluss_clipped(relid, rho2.ring, rho2.weno.ref(1));
#ifdef _A2_A2R2_R_
    _zweno_minus_clipped(relid, rho.ring, rho.weno.ref(0));
    _zweno_pluss_clipped(relid, rho.ring, rho.weno.ref(1));
#endif
    _zweno_minus_clipped(relid, u.ring, u.weno.ref(0));
    _zweno_pluss_clipped(relid, u.ring, u.weno.ref(1));
    _zweno_minus_clipped(relid, v.ring, v.weno.ref(0));
    _zweno_pluss_clipped(relid, v.ring, v.weno.ref(1));
    _zweno_minus_clipped(relid, w.ring, w.weno.ref(0));
    _zweno_pluss_clipped(relid, w.ring, w.weno.ref(1));
    _zweno_minus_clipped(relid, p.ring, p.weno.ref(0));
    _zweno_pluss_clipped(relid, p.ring, p.weno.ref(1));
    _zweno_minus_clipped(relid, A2.ring, A2.weno.ref(0));
    _zweno_pluss_clipped(relid, A2.ring, A2.weno.ref(1));

#ifdef _A2_A2R2_R_
    _convert_reconstruction();
#endif

#ifdef _RECONPCLIP_
   _pressure_clipping_hllc(_BLOCKSIZE_, A2.weno(0), A2.weno(1),p.weno.ref(0), p.weno.ref(1));
#endif

    // characteristic velocities
    _char_vel_einfeldt(_BLOCKSIZE_, rho1.weno(0), rho1.weno(1), rho2.weno(0), rho2.weno(1), w.weno(0), w.weno(1), p.weno(0), p.weno(1), A2.weno(0), A2.weno(1), charvel.ref(0), charvel.ref(1));

    // intermediate wave
    _char_vel_star(_BLOCKSIZE_, rho1.weno(0), rho1.weno(1), rho2.weno(0), rho2.weno(1), w.weno(0), w.weno(1), p.weno(0), p.weno(1), charvel(0), charvel(1), charvel_star.ref(0));

    // RHS of advection equations
    /* TODO: (fabianw; Thu 29 Oct 2015 08:25:01 PM CET) WARNING: According to
     * the ring, A2m = A2.weno(0) and A2p = A2.weno(-1) [the same is true for
     * p.weno()].  Here it does not matter as long as for A2m and pm the
     * corresponding weno() slices are passed, and vice versa for A2p and pp
     * (which is the case here). */
    /* if (evalXtra) _zextraterm_hllc(w.weno(-2), w.weno(-1), w.weno(0), w.weno(1), A2.weno(-1), A2.weno(0), charvel(-2), charvel(-1), charvel_star(-1), charvel(0), charvel(1), charvel_star(0)); */
    if (evalXtra) _zextraterm_hllc(w.weno(-2), w.weno(-1), w.weno(0), w.weno(1), A2.weno(-1), A2.weno(0), p.weno(-1), p.weno(0), charvel(-2), charvel(-1), charvel_star(-1), charvel(0), charvel(1), charvel_star(0));

    // conservative variables
    _hllc_rho( _BLOCKSIZE_, rho1.weno(0), rho1.weno(1), w.weno(0), w.weno(1), charvel(0), charvel(1), charvel_star(0), rho1.flux.ref());
    _hllc_rho( _BLOCKSIZE_, rho2.weno(0), rho2.weno(1), w.weno(0), w.weno(1), charvel(0), charvel(1), charvel_star(0), rho2.flux.ref());
    _hllc_vel( _BLOCKSIZE_, rho1.weno(0), rho1.weno(1), rho2.weno(0), rho2.weno(1), u.weno(0), u.weno(1), w.weno(0), w.weno(1), charvel(0), charvel(1), charvel_star(0), u.flux.ref());
    _hllc_vel( _BLOCKSIZE_, rho1.weno(0), rho1.weno(1), rho2.weno(0), rho2.weno(1), v.weno(0), v.weno(1), w.weno(0), w.weno(1), charvel(0), charvel(1), charvel_star(0), v.flux.ref());
    _hllc_pvel(_BLOCKSIZE_, rho1.weno(0), rho1.weno(1), rho2.weno(0), rho2.weno(1), w.weno(0), w.weno(1), p.weno(0), p.weno(1), charvel(0), charvel(1), charvel_star(0), w.flux.ref());
    _hllc_e(   _BLOCKSIZE_, rho1.weno(0), rho1.weno(1), rho2.weno(0), rho2.weno(1), w.weno(0), w.weno(1), u.weno(0), u.weno(1), v.weno(0), v.weno(1), p.weno(0), p.weno(1), A2.weno(0), A2.weno(1), charvel(0), charvel(1), charvel_star(0), p.flux.ref());

    // advected quantities
    _hllc_rho(_BLOCKSIZE_, A2.weno(0), A2.weno(1), w.weno(0), w.weno(1), charvel(0), charvel(1), charvel_star(0), A2.flux.ref());
}


void Convection_CPP_HLLC_5eq::_xrhs()
{
    _xdivergence(rho1.flux(), rho1.rhs);
    _xdivergence(rho2.flux(), rho2.rhs);
    _xdivergence(u.flux(), u.rhs);
    _xdivergence(v.flux(), v.rhs);
    _xdivergence(w.flux(), w.rhs);
    _xdivergence(p.flux(), p.rhs);
    _xdivergence(A2.flux(), A2.rhs);
}

void Convection_CPP_HLLC_5eq::_yrhs()
{
    _ydivergence(rho1.flux(), rho1.rhs);
    _ydivergence(rho2.flux(), rho2.rhs);
    _ydivergence(u.flux(), u.rhs);
    _ydivergence(v.flux(), v.rhs);
    _ydivergence(w.flux(), w.rhs);
    _ydivergence(p.flux(), p.rhs);
    _ydivergence(A2.flux(), A2.rhs);
}

void Convection_CPP_HLLC_5eq::_zrhs()
{
    _zdivergence(rho1.flux(-1), rho1.flux(0), rho1.rhs);
    _zdivergence(rho2.flux(-1), rho2.flux(0), rho2.rhs);
    _zdivergence(u.flux(-1), u.flux(0), u.rhs);
    _zdivergence(v.flux(-1), v.flux(0), v.rhs);
    _zdivergence(w.flux(-1), w.flux(0), w.rhs);
    _zdivergence(p.flux(-1), p.flux(0), p.rhs);
    _zdivergence(A2.flux(-1), A2.flux(0), A2.rhs);
}


void Convection_CPP_HLLC_5eq::_copyback(Real * const gptfirst, const int gptfloats, const int rowgpts)
{
    const Real factor2 = ((Real)1.)/6;

    for(int iy=0; iy<OutputSOA::NY; iy++)
        for(int ix=0; ix<OutputSOA::NX; ix++)
        {
            AssumedType& rhs = *(AssumedType*)(gptfirst + gptfloats*(ix + iy*rowgpts));

            assert(!std::isnan(rho1.rhs(ix, iy)));
            assert(!std::isnan(rho2.rhs(ix, iy)));
            assert(!std::isnan(u.rhs(ix, iy)));
            assert(!std::isnan(v.rhs(ix, iy)));
            assert(!std::isnan(w.rhs(ix, iy)));
            assert(!std::isnan(p.rhs(ix, iy)));
            assert(!std::isnan(A2.rhs(ix, iy)));
            assert(!std::isnan(sumA2(ix, iy)));

            assert(!std::isinf(rho1.rhs(ix, iy)));
            assert(!std::isinf(rho2.rhs(ix, iy)));
            assert(!std::isinf(u.rhs(ix, iy)));
            assert(!std::isinf(v.rhs(ix, iy)));
            assert(!std::isinf(w.rhs(ix, iy)));
            assert(!std::isinf(p.rhs(ix, iy)));
            assert(!std::isinf(A2.rhs(ix, iy)));
            assert(!std::isinf(sumA2(ix, iy)));
#ifndef _NOK_
            assert(!std::isnan(sumK(ix, iy)));
            assert(!std::isinf(sumK(ix, iy)));
#endif

            rhs.r1 = a*rhs.r1 - dtinvh*rho1.rhs(ix, iy);
            rhs.r2 = a*rhs.r2 - dtinvh*rho2.rhs(ix, iy);
            rhs.u  = a*rhs.u -  dtinvh*u.rhs(ix, iy);
            rhs.v  = a*rhs.v -  dtinvh*v.rhs(ix, iy);
            rhs.w  = a*rhs.w -  dtinvh*w.rhs(ix, iy);
            rhs.E  = a*rhs.E -  dtinvh*p.rhs(ix, iy);
#ifdef _NOK_
            rhs.A2 = a*rhs.A2 - dtinvh*(A2.rhs(ix, iy) - divu(ix,iy)*sumA2(ix,iy)*factor2);
#else
            rhs.A2 = a*rhs.A2 - dtinvh*(A2.rhs(ix, iy) - divu(ix,iy)*sumA2(ix,iy)*factor2 - divu(ix,iy)*sumK(ix,iy)*factor2);
#endif
        }
}
