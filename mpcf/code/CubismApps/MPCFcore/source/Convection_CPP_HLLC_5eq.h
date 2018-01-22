/*
 *  Convection_CPP_HLLC_5eq.h
 *  MPCFcore
 *
 *  Created by Fabian Wermelinger 10/21/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef CONVECTION_CPP_HLLC_5EQ_H_VQUUSOF5
#define CONVECTION_CPP_HLLC_5EQ_H_VQUUSOF5

#include "Convection_CPP.h"

class Convection_CPP_HLLC_5eq : public Convection_CPP
{
public:
    Convection_CPP_HLLC_5eq(const Real a, const Real dtinvh, const Real g1=1.4, const Real g2=1.4, const Real pc1=0.0, const Real pc2=0.0) :
        Convection_CPP(a, dtinvh), _g1(g1), _g2(g2), _pc1(pc1), _pc2(pc2), _g1m1Inv(1.0/(g1-1.0)), _g2m1Inv(1.0/(g2-1.0))
#ifdef _QUANTIFY_DELTAP_
        , _pcorrMax(0), _pcorrCount(0)
#endif /* _QUANTIFY_DELTAP_ */
        {}

#ifdef _QUANTIFY_DELTAP_
    inline Real getMaxPCorrection() const { return _pcorrMax; }
    inline size_t getNumberOfPCorrections() const { return _pcorrCount; }
    inline void resetPCorrection() { _pcorrMax=0; _pcorrCount=0; }
#endif /* _QUANTIFY_DELTAP_ */

protected:

    struct AssumedType { Real r1, r2, u, v, w, E, A2, dummy; };

    const Real _g1, _g2;
    const Real _pc1, _pc2;
    const Real _g1m1Inv, _g2m1Inv; // inverse of constant (gamma_k - 1) for k = {1,2}

#ifdef _QUANTIFY_DELTAP_
    Real _pcorrMax;
    size_t _pcorrCount;
#endif /* _QUANTIFY_DELTAP_ */


    /* TODO: (fabianw; Wed 21 Oct 2015 07:13:55 PM CEST) we just duplicate
     * these attributes here for clarity.  This need refactoring later on. */
    WorkingSet<2> rho1, rho2, u, v;
    WorkingSet<4> w, A2, p;
    OutputSOA sumA2, sumK, divu;
    /* TODO: (fabianw; Wed 21 Oct 2015 07:15:10 PM CEST) stop here */

    // need additional ring with two slices for the velocities of the
    // contact wave
    RingTempSOA charvel_star;

    virtual void _next()
    {
        rho1.ring.next(); rho2.ring.next();
#ifdef _A2_A2R2_R_
        rho.ring.next();
#endif
        u.ring.next(); v.ring.next(); w.ring.next(); p.ring.next(); A2.ring.next();
    }

    virtual void _flux_next()
    {
        rho1.flux.next(); rho2.flux.next(); u.flux.next(); v.flux.next(); w.flux.next(); p.flux.next(); A2.flux.next();
        w.weno.next(); A2.weno.next(); w.weno.next(); A2.weno.next();
        p.weno.next(); p.weno.next();
        charvel.next(); charvel.next();
        charvel_star.next();
    }

#ifdef _RECONPCLIP_
    // helpers
    void _pressure_clipping_hllc(const int SizeX,
            const TempSOA& A2m, const TempSOA& A2p,
            TempSOA& pm, TempSOA& pp);
#endif /* _RECONPCLIP_ */

#ifdef _A2_A2R2_R_
    inline void _convert_reconstruction()
    {
        rho1.weno.ref(0) = rho.weno.ref(0);
        rho1.weno.ref(0) -= rho2.weno.ref(0);
        rho1.weno.ref(1) = rho.weno.ref(1);
        rho1.weno.ref(1) -= rho2.weno.ref(1);
    }
#endif /* _A2_A2R2_R_ */

    // new methods
    void _xextraterm_hllc( const TempSOA& um, const TempSOA& up,
            const TempSOA& A2m, const TempSOA& A2p,
            const TempSOA& pm, const TempSOA& pp,
            const TempSOA& sm, const TempSOA& sp, const TempSOA& ss);

    void _yextraterm_hllc( const TempSOA& vm, const TempSOA& vp,
            const TempSOA& A2m, const TempSOA& A2p,
            const TempSOA& pm, const TempSOA& pp,
            const TempSOA& sm, const TempSOA& sp, const TempSOA& ss);

    void _zextraterm_hllc( const TempSOA& wm0, const TempSOA& wp0,
            const TempSOA& wm1, const TempSOA& wp1,
            const TempSOA& A2m,  const TempSOA& A2p,
            const TempSOA& pm, const TempSOA& pp,
            const TempSOA& sm0, const TempSOA& sp0, const TempSOA& ss0,
            const TempSOA& sm1, const TempSOA& sp1, const TempSOA& ss1);

    // ====================================================================
    void _char_vel_einfeldt( const int SizeX,
            const TempSOA& r1m, const TempSOA& r1p,
            const TempSOA& r2m, const TempSOA& r2p,
            const TempSOA& vm, const TempSOA& vp,
            const TempSOA& pm, const TempSOA& pp,
            const TempSOA& A2m, const TempSOA& A2p,
            TempSOA& out_minus, TempSOA& out_plus);

    void _char_vel_star( const int SizeX,
            const TempSOA& r1m, const TempSOA& r1p,
            const TempSOA& r2m, const TempSOA& r2p,
            const TempSOA& vm, const TempSOA& vp,
            const TempSOA& pm, const TempSOA& pp,
            const TempSOA& sm, const TempSOA& sp,
            TempSOA& out_star);
    // ====================================================================

    void _hllc_rho( const int SizeX,
            const TempSOA& rm, const TempSOA& rp,
            const TempSOA& vm, const TempSOA& vp,
            const TempSOA& sm, const TempSOA& sp, const TempSOA& ss,
            TempSOA& flux_out);

    void _hllc_vel( const int SizeX,
            const TempSOA& r1m, const TempSOA& r1p,
            const TempSOA& r2m, const TempSOA& r2p,
            const TempSOA& vm,  const TempSOA& vp,
            const TempSOA& vdm, const TempSOA& vdp,
            const TempSOA& sm,  const TempSOA& sp,  const TempSOA& ss,
            TempSOA& flux_out);

    void _hllc_pvel( const int SizeX,
            const TempSOA& r1m, const TempSOA& r1p,
            const TempSOA& r2m, const TempSOA& r2p,
            const TempSOA& vm, const TempSOA& vp,
            const TempSOA& pm, const TempSOA& pp,
            const TempSOA& sm, const TempSOA& sp, const TempSOA& ss,
            TempSOA& flux_out);

    void _hllc_e( const int SizeX,
            const TempSOA& r1m, const TempSOA& r1p,
            const TempSOA& r2m, const TempSOA& r2p,
            const TempSOA& vdm, const TempSOA& vdp,
            const TempSOA& v1m, const TempSOA& v1p,
            const TempSOA& v2m, const TempSOA& v2p,
            const TempSOA& pm,  const TempSOA& pp,
            const TempSOA& A2m,  const TempSOA& A2p,
            const TempSOA& sm,  const TempSOA& sp,  const TempSOA& ss,
            TempSOA& flux_out);

    // inherited and modified in this class
    virtual void _convert(const Real * const gptfirst, const int gptfloats, const int rowgpts, const int islice);

    virtual void _xflux(const int relsliceid);
    virtual void _yflux(const int relsliceid);
    virtual void _zflux(const int relsliceid, const bool evalXtra);

    virtual void _xrhs();
    virtual void _yrhs();
    virtual void _zrhs();

    virtual void _copyback(Real * const gptfirst, const int gptfloats, const int rowgpts);
};

#endif /* CONVECTION_CPP_HLLC_5EQ_H_VQUUSOF5 */
