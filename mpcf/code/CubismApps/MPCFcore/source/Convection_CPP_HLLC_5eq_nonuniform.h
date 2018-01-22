/*
 *  Convection_CPP_HLLC_5eq_nonuniform.h
 *  MPCFcore
 *
 *  Created by Fabian Wermelinger 05/09/2017
 *  Copyright 2017 ETH Zurich. All rights reserved.
 *
 */
#ifndef CONVECTION_CPP_HLLC_5EQ_NONUNIFORM_H_UX1KPZYB
#define CONVECTION_CPP_HLLC_5EQ_NONUNIFORM_H_UX1KPZYB

#include "WenoCoefficients.h"
#include "Convection_CPP_HLLC_5eq.h"

class Convection_CPP_HLLC_5eq_nonuniform : public Convection_CPP_HLLC_5eq
{
public:
    Convection_CPP_HLLC_5eq_nonuniform(const Real a, const Real dt,
            const Real g1=1.4, const Real g2=1.4,
            const Real pc1=0.0, const Real pc2=0.0) :
        Convection_CPP_HLLC_5eq(a, 0.0, g1,g2,pc1,pc2), // explicitly set dtinvh = 0!
        m_dt(dt)
        {}

    void compute_nonuniform(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
            Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts,
            const Coefficients_t& xcoeffs_minus, const Coefficients_t& xcoeffs_plus,
            const Coefficients_t& ycoeffs_minus, const Coefficients_t& ycoeffs_plus,
            const Coefficients_t& zcoeffs_minus, const Coefficients_t& zcoeffs_plus,
            const Real* const invh_x, const Real* const invh_y, const Real* const invh_z);

protected:
    const Real m_dt;

    void _xweno_minus(const InputSOA& in, TempSOA& out, const Coefficients_t& coeffs);
    void _xweno_pluss(const InputSOA& in, TempSOA& out, const Coefficients_t& coeffs);
    void _yweno_minus(const InputSOA& in, TempSOA& out, const Coefficients_t& coeffs);
    void _yweno_pluss(const InputSOA& in, TempSOA& out, const Coefficients_t& coeffs);
    void _zweno_minus(const int relid, const RingInputSOA& in, TempSOA& out, const int sliceid, const Coefficients_t& coeffs);
    void _zweno_pluss(const int relid, const RingInputSOA& in, TempSOA& out, const int sliceid, const Coefficients_t& coeffs);

    void _xextraterm_hllc_nonuniform(const TempSOA& um, const TempSOA& up,
            const TempSOA& A2m, const TempSOA& A2p,
            const TempSOA& pm, const TempSOA& pp,
            const TempSOA& sm, const TempSOA& sp, const TempSOA& ss,
            const Real* const invh);

    void _yextraterm_hllc_nonuniform(const TempSOA& vm, const TempSOA& vp,
            const TempSOA& A2m, const TempSOA& A2p,
            const TempSOA& pm, const TempSOA& pp,
            const TempSOA& sm, const TempSOA& sp, const TempSOA& ss,
            const Real* const invh);

    void _zextraterm_hllc_nonuniform(const TempSOA& wm0, const TempSOA& wp0,
            const TempSOA& wm1, const TempSOA& wp1,
            const TempSOA& A2m,  const TempSOA& A2p,
            const TempSOA& pm, const TempSOA& pp,
            const TempSOA& sm0, const TempSOA& sp0, const TempSOA& ss0,
            const TempSOA& sm1, const TempSOA& sp1, const TempSOA& ss1,
            const Real invh);

    virtual void _xdivergence_nonuniform(const TempSOA& flux, OutputSOA& rhs, const Real* const invh);
    virtual void _ydivergence_nonuniform(const TempSOA& flux, OutputSOA& rhs, const Real* const invh);
    virtual void _zdivergence_nonuniform(const TempSOA& fback, const TempSOA& fforward, OutputSOA& rhs, const Real invh);

    virtual void _xflux_nonuniform(const int relsliceid, const Coefficients_t& coeffs_minus, const Coefficients_t& coeffs_plus, const Real* const invh);
    virtual void _yflux_nonuniform(const int relsliceid, const Coefficients_t& coeffs_minus, const Coefficients_t& coeffs_plus, const Real* const invh);
    virtual void _zflux_nonuniform(const int relsliceid, const bool evalXtra, const int sliceid, const Coefficients_t& coeffs_minus, const Coefficients_t& coeffs_plus, const Real invh);

    virtual void _xrhs_nonuniform(const Real* const invh);
    virtual void _yrhs_nonuniform(const Real* const invh);
    virtual void _zrhs_nonuniform(const Real invh);

    virtual void _copyback(Real * const gptfirst, const int gptfloats, const int rowgpts);
};

#endif /* CONVECTION_CPP_HLLC_5EQ_NONUNIFORM_H_UX1KPZYB */
