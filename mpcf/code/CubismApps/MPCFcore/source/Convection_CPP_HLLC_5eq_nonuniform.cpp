/*
 *  Convection_CPP_HLLC_5eq_nonuniform.cpp
 *  MPCFcore
 *
 *  Created by Fabian Wermelinger 05/09/2017
 *  Copyright 2017 ETH Zurich. All rights reserved.
 *
 */
#include <cassert>
#include <cmath>
#include <cstdio>
#include "Convection_CPP_HLLC_5eq_nonuniform.h"

using namespace std;

// NOTE:
// oefficients minus correspond to the following sketch:
// |  a  |  b  |  c -|+ d  |  e  |  f  |
//
// Therefore, coefficients evaluated in cell c at -| correspond to location
// x_{i+1/2} in that cell and vice versa in cell d.
void Convection_CPP_HLLC_5eq_nonuniform::compute_nonuniform(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
        Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts,
        const Coefficients_t& xcoeffs_minus, const Coefficients_t& xcoeffs_plus,
        const Coefficients_t& ycoeffs_minus, const Coefficients_t& ycoeffs_plus,
        const Coefficients_t& zcoeffs_minus, const Coefficients_t& zcoeffs_plus,
        const Real* const invh_x, const Real* const invh_y, const Real* const invh_z)
{
    for(int islice=0; islice<5; islice++)
    {
        _convert(srcfirst+islice*srcfloats*slicesrcs, srcfloats, rowsrcs, islice);
        _next();
    }

    _convert(srcfirst + 5*srcfloats*slicesrcs, srcfloats, rowsrcs, 5);

    _zflux_nonuniform(-2, false, 0, zcoeffs_minus, zcoeffs_plus, 0.0);  // first pass -> illegal to evaluate _zextraterm(); inverse h_z is not defined
    _flux_next();

    for(int islice=0; islice<_BLOCKSIZE_; islice++)
    {
        _xflux_nonuniform(-2, xcoeffs_minus, xcoeffs_plus, invh_x);
        _xrhs_nonuniform(invh_x);

        _yflux_nonuniform(-2, ycoeffs_minus, ycoeffs_plus, invh_y);
        _yrhs_nonuniform(invh_y);

        _next();
        _convert(srcfirst + (islice+6)*srcfloats*slicesrcs, srcfloats, rowsrcs, islice+6);

        const Real ihz = invh_z[islice];
        _zflux_nonuniform(-2, true, islice+1, zcoeffs_minus, zcoeffs_plus, ihz); // information for previous slice is computed now -> ok to evaluate _zextraterm()
        _zrhs_nonuniform(ihz);

        _copyback(dstfirst + islice*dstfloats*slicedsts, dstfloats, rowdsts);
        _flux_next();
    }
}


///////////////////////////////////////////////////////////////////////////////
// Static methods
///////////////////////////////////////////////////////////////////////////////
// |  a  |  b  |  c -|+ d  |  e  |  f  |
inline Real weno_minus_nonuniform(const Real a, const Real b, const Real c, const Real d, const Real e,
        const Real b00, const Real b01, const Real b02,
        const Real b10, const Real b11, const Real b12,
        const Real b20, const Real b21, const Real b22,
        const Real c00, const Real c01, const Real c02,
        const Real c10, const Real c11, const Real c12,
        const Real c20, const Real c21, const Real c22,
        const Real d0,  const Real d1,  const Real d2)
{
#ifdef _WENO3_
    // compute nonlinear weights
    ///////////////////////////////////////////////////////////////////////////
    // smoothness indicators
    const Real IS0 = b00*c*c + b01*c*d + b02*d*d;
    const Real IS1 = b10*b*b + b11*b*c + b12*c*c;
    assert(IS0>=0.);
    assert(IS1>=0.);

    const Real IS0_den = IS0 + (Real)WENOEPS;
    const Real IS1_den = IS1 + (Real)WENOEPS;

    // scaled ideal weights
    const Real alpha0 = d0 / (IS0_den*IS0_den);
    const Real alpha1 = d1 / (IS1_den*IS1_den);
    const Real inv_alpha = (Real)1.0 / (alpha0 + alpha1);
    assert(alpha0>=0. && !std::isnan(alpha0) && !std::isinf(alpha0));
    assert(alpha1>=0. && !std::isnan(alpha1) && !std::isinf(alpha1));
    assert(inv_alpha>=0. && !std::isnan(inv_alpha) && !std::isinf(inv_alpha));

    // nonlinear weights
    const Real omega0 = alpha0 * inv_alpha;
    // const Real omega1 = alpha1 * inv_alpha;
    const Real omega1 = (Real)1.0 - omega0;
    assert(omega0>=0. && !std::isnan(omega0) && !std::isinf(omega0));
    assert(omega1>=0. && !std::isnan(omega1) && !std::isinf(omega1));

    // compute polynomials
    ///////////////////////////////////////////////////////////////////////////
    const Real p0 = c00*c + c01*d;
    const Real p1 = c10*b + c11*c;

    return omega0*p0 + omega1*p1;
#else /* WENO5 */
    // compute nonlinear weights
    ///////////////////////////////////////////////////////////////////////////
    // smoothness indicators
    const Real IS0 = b00*(d-c)*(d-c) + b01*(d-c)*(e-d) + b02*(e-d)*(e-d);
    const Real IS1 = b10*(c-b)*(d-c) + b11*(c-b)*(c-b) + b12*(d-c)*(d-c);
    const Real IS2 = b20*(c-b)*(c-b) + b21*(b-a)*(c-b) + b22*(b-a)*(b-a);
    assert(IS0>=0.);
    assert(IS1>=0.);
    assert(IS2>=0.);

    const Real IS0_den = IS0 + (Real)WENOEPS;
    const Real IS1_den = IS1 + (Real)WENOEPS;
    const Real IS2_den = IS2 + (Real)WENOEPS;

    // scaled ideal weights
    const Real alpha0 = d0 / (IS0_den*IS0_den);
    const Real alpha1 = d1 / (IS1_den*IS1_den);
    const Real alpha2 = d2 / (IS2_den*IS2_den);
    const Real inv_alpha = (Real)1.0 / (alpha0 + alpha1 + alpha2);
    assert(alpha0>=0. && !std::isnan(alpha0) && !std::isinf(alpha0));
    assert(alpha1>=0. && !std::isnan(alpha1) && !std::isinf(alpha1));
    assert(alpha2>=0. && !std::isnan(alpha2) && !std::isinf(alpha2));
    assert(inv_alpha>=0. && !std::isnan(inv_alpha) && !std::isinf(inv_alpha));

    // nonlinear weights
    const Real omega0 = alpha0 * inv_alpha;
    const Real omega1 = alpha1 * inv_alpha;
    // const Real omega2 = alpha2 * inv_alpha;
    const Real omega2 = (Real)1.0 - omega0 - omega1;
    assert(omega0>=0. && !std::isnan(omega0) && !std::isinf(omega0));
    assert(omega1>=0. && !std::isnan(omega1) && !std::isinf(omega1));
    assert(omega2>=0. && !std::isnan(omega2) && !std::isinf(omega2));

    // compute polynomials (based on Coralic)
    ///////////////////////////////////////////////////////////////////////////
    const Real p0 = c + c01*(d-c) + c02*(e-d);
    const Real p1 = c + c11*(c-b) + c12*(d-c);
    const Real p2 = c + c21*(b-a) + c22*(c-b);

    return omega0*p0 + omega1*p1 + omega2*p2;
#endif /* _WENO3_ */
}


// |  a  |  b  |  c -|+ d  |  e  |  f  |
inline Real weno_plus_nonuniform(const Real b, const Real c, const Real d, const Real e, const Real f,
        const Real b00, const Real b01, const Real b02,
        const Real b10, const Real b11, const Real b12,
        const Real b20, const Real b21, const Real b22,
        const Real c00, const Real c01, const Real c02,
        const Real c10, const Real c11, const Real c12,
        const Real c20, const Real c21, const Real c22,
        const Real d0,  const Real d1,  const Real d2)
{
#ifdef _WENO3_
    // compute nonlinear weights
    ///////////////////////////////////////////////////////////////////////////
    // smoothness indicators
    const Real IS0 = b00*d*d + b01*d*e + b02*e*e;
    const Real IS1 = b10*c*c + b11*c*d + b12*d*d;
    assert(IS0>=0.);
    assert(IS1>=0.);

    const Real IS0_den = IS0 + (Real)WENOEPS;
    const Real IS1_den = IS1 + (Real)WENOEPS;

    // scaled ideal weights
    const Real alpha0 = d0 / (IS0_den*IS0_den);
    const Real alpha1 = d1 / (IS1_den*IS1_den);
    const Real inv_alpha = (Real)1.0 / (alpha0 + alpha1);
    assert(alpha0>=0. && !std::isnan(alpha0) && !std::isinf(alpha0));
    assert(alpha1>=0. && !std::isnan(alpha1) && !std::isinf(alpha1));
    assert(inv_alpha>=0. && !std::isnan(inv_alpha) && !std::isinf(inv_alpha));

    // nonlinear weights
    const Real omega0 = alpha0 * inv_alpha;
    // const Real omega1 = alpha1 * inv_alpha;
    const Real omega1 = (Real)1.0 - omega0;
    assert(omega0>=0. && !std::isnan(omega0) && !std::isinf(omega0));
    assert(omega1>=0. && !std::isnan(omega1) && !std::isinf(omega1));

    // compute polynomials
    ///////////////////////////////////////////////////////////////////////////
    const Real p0 = c00*d + c01*e;
    const Real p1 = c10*c + c11*d;

    return omega0*p0 + omega1*p1;
#else /* WENO5 */
    // compute nonlinear weights
    ///////////////////////////////////////////////////////////////////////////
    // smoothness indicators
    const Real IS0 = b00*(e-d)*(e-d) + b01*(e-d)*(f-e) + b02*(f-e)*(f-e);
    const Real IS1 = b10*(d-c)*(e-d) + b11*(d-c)*(d-c) + b12*(e-d)*(e-d);
    const Real IS2 = b20*(d-c)*(d-c) + b21*(c-b)*(d-c) + b22*(c-b)*(c-b);
    assert(IS0>=0.);
    assert(IS1>=0.);
    assert(IS2>=0.);

    const Real IS0_den = IS0 + (Real)WENOEPS;
    const Real IS1_den = IS1 + (Real)WENOEPS;
    const Real IS2_den = IS2 + (Real)WENOEPS;

    // scaled ideal weights
    const Real alpha0 = d0 / (IS0_den*IS0_den);
    const Real alpha1 = d1 / (IS1_den*IS1_den);
    const Real alpha2 = d2 / (IS2_den*IS2_den);
    const Real inv_alpha = (Real)1.0 / (alpha0 + alpha1 + alpha2);
    assert(alpha0>=0. && !std::isnan(alpha0) && !std::isinf(alpha0));
    assert(alpha1>=0. && !std::isnan(alpha1) && !std::isinf(alpha1));
    assert(alpha2>=0. && !std::isnan(alpha2) && !std::isinf(alpha2));
    assert(inv_alpha>=0. && !std::isnan(inv_alpha) && !std::isinf(inv_alpha));

    // nonlinear weights
    const Real omega0 = alpha0 * inv_alpha;
    const Real omega1 = alpha1 * inv_alpha;
    // const Real omega2 = alpha2 * inv_alpha;
    const Real omega2 = (Real)1.0 - omega0 - omega1;
    assert(omega0>=0. && !std::isnan(omega0) && !std::isinf(omega0));
    assert(omega1>=0. && !std::isnan(omega1) && !std::isinf(omega1));
    assert(omega2>=0. && !std::isnan(omega2) && !std::isinf(omega2));

    // compute polynomials (based on Coralic)
    ///////////////////////////////////////////////////////////////////////////
    const Real p0 = d + c01*(e-d) + c02*(f-e);
    const Real p1 = d + c11*(d-c) + c12*(e-d);
    const Real p2 = d + c21*(c-b) + c22*(d-c);

    return omega0*p0 + omega1*p1 + omega2*p2;
#endif /* _WENO3_ */
}


inline Real weno_minus_nonuniform_clipped(const Real a, const Real b, const Real c, const Real d, const Real e,
        const Real b00, const Real b01, const Real b02,
        const Real b10, const Real b11, const Real b12,
        const Real b20, const Real b21, const Real b22,
        const Real c00, const Real c01, const Real c02,
        const Real c10, const Real c11, const Real c12,
        const Real c20, const Real c21, const Real c22,
        const Real d0,  const Real d1,  const Real d2)
{
       const Real retval = weno_minus_nonuniform(a,b,c,d,e,
               b00, b01, b02, b10, b11, b12, b20, b21, b22,
               c00, c01, c02, c10, c11, c12, c20, c21, c22,
               d0, d1, d2);

#ifndef _ADCLIP_

// #ifdef _WENO3_
//        const Real min_in = min(min(b,c),d);
//        const Real max_in = max(max(b,c),d);
// #else
       const Real min_in = min(min(a,b),min(min(c,d),e));
       const Real max_in = max(max(a,b),max(max(c,d),e));
// #endif

       return min(max((Real)retval, min_in), max_in);

#else
       const Real di = d - 2.0*c + b;
       const Real dip1 = e - 2.0*d + c;
       const Real dim1 = c - 2.0*b + a;

       const Real diff = 4.0*di - dip1;
       const Real diffp = 4.0*dip1-di;
       const Real dm4ip = 0.5*(copysign(1.0,diff) + copysign(1.0,diffp)) * abs(0.25 * (copysign(1.0,diff)+copysign(1.0,di)) * (copysign(1.0,diff)+copysign(1.0,dip1))) * min(min(min(abs(diff),abs(diffp)),abs(di)),abs(dip1));

       const Real diffl = 4.0*dim1 - di;
       const Real diffm = 4.0*di-dim1;
       const Real dm4im = 0.5*(copysign(1.0,diffl) + copysign(1.0,diffm)) * abs(0.25 * (copysign(1.0,diffl)+copysign(1.0,di)) * (copysign(1.0,diffl)+copysign(1.0,dim1))) * min(min(min(abs(diffl),abs(diffm)),abs(di)),abs(dim1));

       const Real u_UL = c + 2.0 * (c-b);
       const Real u_MD = 0.5 * (c+d) - 0.5 * dm4ip;
       const Real u_LC = c + 0.5 * (c-b) + 4.0/3.0 * dm4im;

       const Real u_minus_min = max(min(min(c,d),u_MD),min(min(c,u_UL),u_LC));
       const Real u_minus_max = min(max(max(c,d),u_MD),max(max(c,u_UL),u_LC));

       const Real diff_min = u_minus_min - retval;
       const Real diff_max = u_minus_max - retval;

       const Real val = (retval + 0.5*(copysign(1.0,diff_min) + copysign(1.0,diff_max)) * min(abs(diff_min), abs(diff_max)));
       //return (retval + 0.5*(diff_min/abs(diff_min) + diff_max/abs(diff_max)) * min(abs(diff_min), abs(diff_max)));

       return val;
#endif
}

inline Real weno_plus_nonuniform_clipped(const Real b, const Real c, const Real d, const Real e, const Real f,
        const Real b00, const Real b01, const Real b02,
        const Real b10, const Real b11, const Real b12,
        const Real b20, const Real b21, const Real b22,
        const Real c00, const Real c01, const Real c02,
        const Real c10, const Real c11, const Real c12,
        const Real c20, const Real c21, const Real c22,
        const Real d0,  const Real d1,  const Real d2)
{
       const Real retval = weno_plus_nonuniform(b,c,d,e,f,
               b00, b01, b02, b10, b11, b12, b20, b21, b22,
               c00, c01, c02, c10, c11, c12, c20, c21, c22,
               d0, d1, d2);

#ifndef _ADCLIP_

// #ifdef _WENO3_
//        const Real min_in = min(min(c,d),e);
//        const Real max_in = max(max(c,d),e);
// #else
       const Real min_in = min(min(b,c),min(min(d,e),f));
       const Real max_in = max(max(b,c),max(max(d,e),f));
// #endif

       return min(max((Real)retval, min_in), max_in);
#else
       const Real di = e - 2.0*d + c;
       const Real dip1 = f - 2.0*e + d;
       const Real dim1 = d - 2.0*c + b;

       const Real diff = 4.0*di - dim1;
       const Real diffm = 4.0*dim1-di;
       const Real dm4im = 0.5*(copysign(1.0,diff) + copysign(1.0,diffm)) * abs(0.25 * (copysign(1.0,diff)+copysign(1.0,di)) * (copysign(1.0,diff)+copysign(1.0,dim1))) * min(min(min(abs(diff),abs(diffm)),abs(di)),abs(dim1));

       const Real diffu = 4.0*dip1 - di;
       const Real diffp = 4.0*di-dip1;
       const Real dm4ip = 0.5*(copysign(1.0,diffu) + copysign(1.0,diffp)) * abs(0.25 * (copysign(1.0,diffu)+copysign(1.0,di)) * (copysign(1.0,diffu)+copysign(1.0,dip1))) * min(min(min(abs(diffu),abs(diffp)),abs(di)),abs(dip1));

       const Real u_UL = d + 2.0 * (d-e);
       const Real u_MD = 0.5 * (c+d) - 0.5 * dm4im;
       const Real u_LC = d + 0.5 * (d-e) + 4.0/3.0 * dm4ip;

       const Real u_minus_min = max(min(min(c,d),u_MD),min(min(d,u_UL),u_LC));
       const Real u_minus_max = min(max(max(c,d),u_MD),max(max(d,u_UL),u_LC));

       const Real diff_min = u_minus_min - retval;
       const Real diff_max = u_minus_max - retval;

       const Real val = (retval + 0.5*(copysign(1.0,diff_min) + copysign(1.0,diff_max)) * min(abs(diff_min), abs(diff_max)));
       //return (retval + 0.5*(diff_min/abs(diff_min) + diff_max/abs(diff_max)) * min(diff_min, diff_max));

       return val;
#endif
}

typedef Real (*_weno_t)(const Real a, const Real b, const Real c, const Real d, const Real e,
        const Real b00, const Real b01, const Real b02,
        const Real b10, const Real b11, const Real b12,
        const Real b20, const Real b21, const Real b22,
        const Real c00, const Real c01, const Real c02,
        const Real c10, const Real c11, const Real c12,
        const Real c20, const Real c21, const Real c22,
        const Real d0,  const Real d1,  const Real d2);

#ifdef _WENO_UNCLIP_
_weno_t weno_minus = &weno_minus_nonuniform;
_weno_t weno_plus  = &weno_plus_nonuniform;
#else
_weno_t weno_minus = &weno_minus_nonuniform_clipped;
_weno_t weno_plus  = &weno_plus_nonuniform_clipped;
#endif /* _WENO_UNCLIP_ */
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// Methods
///////////////////////////////////////////////////////////////////////////////
void Convection_CPP_HLLC_5eq_nonuniform::_xweno_minus(const InputSOA& _in, TempSOA& _out, const Coefficients_t& coeffs)
{
    const Real* const c00 = &coeffs.c[0][0];
    const Real* const c01 = &coeffs.c[1][0];
    const Real* const c02 = &coeffs.c[2][0];
    const Real* const c10 = &coeffs.c[3][0];
    const Real* const c11 = &coeffs.c[4][0];
    const Real* const c12 = &coeffs.c[5][0];
    const Real* const c20 = &coeffs.c[6][0];
    const Real* const c21 = &coeffs.c[7][0];
    const Real* const c22 = &coeffs.c[8][0];

    const Real* const d0 = &coeffs.d[0][0];
    const Real* const d1 = &coeffs.d[1][0];
    const Real* const d2 = &coeffs.d[2][0];

    const Real* const b00 = &coeffs.b[0][0];
    const Real* const b01 = &coeffs.b[1][0];
    const Real* const b02 = &coeffs.b[2][0];
    const Real* const b10 = &coeffs.b[3][0];
    const Real* const b11 = &coeffs.b[4][0];
    const Real* const b12 = &coeffs.b[5][0];
    const Real* const b20 = &coeffs.b[6][0];
    const Real* const b21 = &coeffs.b[7][0];
    const Real* const b22 = &coeffs.b[8][0];

    for(int iy=0; iy<TempSOA::NY; iy++)
    {
        const Real * const in = _in.ptr(-3, iy);
        Real * const out = & _out.ref(0, iy);

        for(int ix=0; ix<TempSOA::NX; ix++)
            out[ix] = weno_minus(in[ix+0], in[ix+1], in[ix+2], in[ix+3], in[ix+4],
                    b00[ix], b01[ix], b02[ix], b10[ix], b11[ix], b12[ix], b20[ix], b21[ix], b22[ix],
                    c00[ix], c01[ix], c02[ix], c10[ix], c11[ix], c12[ix], c20[ix], c21[ix], c22[ix],
                    d0[ix], d1[ix], d2[ix]);
    }
}

void Convection_CPP_HLLC_5eq_nonuniform::_xweno_pluss(const InputSOA& _in, TempSOA& _out, const Coefficients_t& coeffs)
{
    const Real* const c00 = &coeffs.c[0][0];
    const Real* const c01 = &coeffs.c[1][0];
    const Real* const c02 = &coeffs.c[2][0];
    const Real* const c10 = &coeffs.c[3][0];
    const Real* const c11 = &coeffs.c[4][0];
    const Real* const c12 = &coeffs.c[5][0];
    const Real* const c20 = &coeffs.c[6][0];
    const Real* const c21 = &coeffs.c[7][0];
    const Real* const c22 = &coeffs.c[8][0];

    const Real* const d0 = &coeffs.d[0][0];
    const Real* const d1 = &coeffs.d[1][0];
    const Real* const d2 = &coeffs.d[2][0];

    const Real* const b00 = &coeffs.b[0][0];
    const Real* const b01 = &coeffs.b[1][0];
    const Real* const b02 = &coeffs.b[2][0];
    const Real* const b10 = &coeffs.b[3][0];
    const Real* const b11 = &coeffs.b[4][0];
    const Real* const b12 = &coeffs.b[5][0];
    const Real* const b20 = &coeffs.b[6][0];
    const Real* const b21 = &coeffs.b[7][0];
    const Real* const b22 = &coeffs.b[8][0];

    for(int iy=0; iy<TempSOA::NY; iy++)
    {
        const Real * const in = _in.ptr(-2, iy);
        Real * const out = & _out.ref(0, iy);

        for(int ix=0; ix<TempSOA::NX; ix++)
            out[ix] = weno_plus(in[ix+0], in[ix+1], in[ix+2], in[ix+3], in[ix+4],
                    b00[ix], b01[ix], b02[ix], b10[ix], b11[ix], b12[ix], b20[ix], b21[ix], b22[ix],
                    c00[ix], c01[ix], c02[ix], c10[ix], c11[ix], c12[ix], c20[ix], c21[ix], c22[ix],
                    d0[ix], d1[ix], d2[ix]);
    }
}

void Convection_CPP_HLLC_5eq_nonuniform::_yweno_minus(const InputSOA& _in, TempSOA& _out, const Coefficients_t& coeffs)
{
    const Real* const c00 = &coeffs.c[0][0];
    const Real* const c01 = &coeffs.c[1][0];
    const Real* const c02 = &coeffs.c[2][0];
    const Real* const c10 = &coeffs.c[3][0];
    const Real* const c11 = &coeffs.c[4][0];
    const Real* const c12 = &coeffs.c[5][0];
    const Real* const c20 = &coeffs.c[6][0];
    const Real* const c21 = &coeffs.c[7][0];
    const Real* const c22 = &coeffs.c[8][0];

    const Real* const d0 = &coeffs.d[0][0];
    const Real* const d1 = &coeffs.d[1][0];
    const Real* const d2 = &coeffs.d[2][0];

    const Real* const b00 = &coeffs.b[0][0];
    const Real* const b01 = &coeffs.b[1][0];
    const Real* const b02 = &coeffs.b[2][0];
    const Real* const b10 = &coeffs.b[3][0];
    const Real* const b11 = &coeffs.b[4][0];
    const Real* const b12 = &coeffs.b[5][0];
    const Real* const b20 = &coeffs.b[6][0];
    const Real* const b21 = &coeffs.b[7][0];
    const Real* const b22 = &coeffs.b[8][0];

    static const int L = InputSOA::PITCH;
    const Real * const in = _in.ptr(0,-3);
    Real * out = &_out.ref(0,0);

    for(int iy=0; iy<TempSOA::NY; iy++)
    {
        const Real * ptr = &in[iy];

        for(int ix=0; ix<TempSOA::NX; ix++)
            out[ix + iy*TempSOA::PITCH] = weno_minus(ptr[ix*L], ptr[ix*L+L], ptr[ix*L+2*L], ptr[ix*L+3*L], ptr[ix*L+4*L],
                    b00[ix], b01[ix], b02[ix], b10[ix], b11[ix], b12[ix], b20[ix], b21[ix], b22[ix],
                    c00[ix], c01[ix], c02[ix], c10[ix], c11[ix], c12[ix], c20[ix], c21[ix], c22[ix],
                    d0[ix], d1[ix], d2[ix]);
    }
}

void Convection_CPP_HLLC_5eq_nonuniform::_yweno_pluss(const InputSOA& _in, TempSOA& _out, const Coefficients_t& coeffs)
{
    const Real* const c00 = &coeffs.c[0][0];
    const Real* const c01 = &coeffs.c[1][0];
    const Real* const c02 = &coeffs.c[2][0];
    const Real* const c10 = &coeffs.c[3][0];
    const Real* const c11 = &coeffs.c[4][0];
    const Real* const c12 = &coeffs.c[5][0];
    const Real* const c20 = &coeffs.c[6][0];
    const Real* const c21 = &coeffs.c[7][0];
    const Real* const c22 = &coeffs.c[8][0];

    const Real* const d0 = &coeffs.d[0][0];
    const Real* const d1 = &coeffs.d[1][0];
    const Real* const d2 = &coeffs.d[2][0];

    const Real* const b00 = &coeffs.b[0][0];
    const Real* const b01 = &coeffs.b[1][0];
    const Real* const b02 = &coeffs.b[2][0];
    const Real* const b10 = &coeffs.b[3][0];
    const Real* const b11 = &coeffs.b[4][0];
    const Real* const b12 = &coeffs.b[5][0];
    const Real* const b20 = &coeffs.b[6][0];
    const Real* const b21 = &coeffs.b[7][0];
    const Real* const b22 = &coeffs.b[8][0];

    static const int L = InputSOA::PITCH;
    const Real * const in = _in.ptr(0,-2);
    Real * out = &_out.ref(0,0);

    for(int iy=0; iy<TempSOA::NY; iy++)
    {
        const Real * ptr = &in[iy];

        for(int ix=0; ix<TempSOA::NX; ix++)
            out[ix + iy*TempSOA::PITCH] = weno_plus(ptr[ix*L], ptr[ix*L+L], ptr[ix*L+2*L], ptr[ix*L+3*L], ptr[ix*L+4*L],
                    b00[ix], b01[ix], b02[ix], b10[ix], b11[ix], b12[ix], b20[ix], b21[ix], b22[ix],
                    c00[ix], c01[ix], c02[ix], c10[ix], c11[ix], c12[ix], c20[ix], c21[ix], c22[ix],
                    d0[ix], d1[ix], d2[ix]);
    }
}

void Convection_CPP_HLLC_5eq_nonuniform::_zweno_minus(const int r, const RingInputSOA& in, TempSOA& out, const int sliceid, const Coefficients_t& coeffs)
{
    static const int L = InputSOA::PITCH;

    const Real * const a = in(r-3).ptr(0,0);
    const Real * const b = in(r-2).ptr(0,0);
    const Real * const c = in(r-1).ptr(0,0);
    const Real * const d = in(r).ptr(0,0);
    const Real * const e = in(r+1).ptr(0,0);

    Real * const o = &out.ref(0,0);

    const Real c00 = coeffs.c[0][sliceid];
    const Real c01 = coeffs.c[1][sliceid];
    const Real c02 = coeffs.c[2][sliceid];
    const Real c10 = coeffs.c[3][sliceid];
    const Real c11 = coeffs.c[4][sliceid];
    const Real c12 = coeffs.c[5][sliceid];
    const Real c20 = coeffs.c[6][sliceid];
    const Real c21 = coeffs.c[7][sliceid];
    const Real c22 = coeffs.c[8][sliceid];

    const Real d0 = coeffs.d[0][sliceid];
    const Real d1 = coeffs.d[1][sliceid];
    const Real d2 = coeffs.d[2][sliceid];

    const Real b00 = coeffs.b[0][sliceid];
    const Real b01 = coeffs.b[1][sliceid];
    const Real b02 = coeffs.b[2][sliceid];
    const Real b10 = coeffs.b[3][sliceid];
    const Real b11 = coeffs.b[4][sliceid];
    const Real b12 = coeffs.b[5][sliceid];
    const Real b20 = coeffs.b[6][sliceid];
    const Real b21 = coeffs.b[7][sliceid];
    const Real b22 = coeffs.b[8][sliceid];

    for(int iy=0; iy<TempSOA::NY; iy++)
        for(int ix=0; ix<TempSOA::NX-1; ix++)
            o[ix + TempSOA::PITCH*iy] = weno_minus(a[ix+L*iy], b[ix+L*iy], c[ix+L*iy], d[ix+L*iy], e[ix+L*iy],
                    b00, b01, b02, b10, b11, b12, b20, b21, b22,
                    c00, c01, c02, c10, c11, c12, c20, c21, c22,
                    d0, d1, d2);
}

void Convection_CPP_HLLC_5eq_nonuniform::_zweno_pluss(const int r, const RingInputSOA& in, TempSOA& out, const int sliceid, const Coefficients_t& coeffs)
{
    static const int L = InputSOA::PITCH;

    const Real * const a = in(r-2).ptr(0,0);
    const Real * const b = in(r-1).ptr(0,0);
    const Real * const c = in(r).ptr(0,0);
    const Real * const d = in(r+1).ptr(0,0);
    const Real * const e = in(r+2).ptr(0,0);

    Real * const o = &out.ref(0,0);

    const Real c00 = coeffs.c[0][sliceid];
    const Real c01 = coeffs.c[1][sliceid];
    const Real c02 = coeffs.c[2][sliceid];
    const Real c10 = coeffs.c[3][sliceid];
    const Real c11 = coeffs.c[4][sliceid];
    const Real c12 = coeffs.c[5][sliceid];
    const Real c20 = coeffs.c[6][sliceid];
    const Real c21 = coeffs.c[7][sliceid];
    const Real c22 = coeffs.c[8][sliceid];

    const Real d0 = coeffs.d[0][sliceid];
    const Real d1 = coeffs.d[1][sliceid];
    const Real d2 = coeffs.d[2][sliceid];

    const Real b00 = coeffs.b[0][sliceid];
    const Real b01 = coeffs.b[1][sliceid];
    const Real b02 = coeffs.b[2][sliceid];
    const Real b10 = coeffs.b[3][sliceid];
    const Real b11 = coeffs.b[4][sliceid];
    const Real b12 = coeffs.b[5][sliceid];
    const Real b20 = coeffs.b[6][sliceid];
    const Real b21 = coeffs.b[7][sliceid];
    const Real b22 = coeffs.b[8][sliceid];

    for(int iy=0; iy<TempSOA::NY; iy++)
        for(int ix=0; ix<TempSOA::NX-1; ix++)
            o[ix + TempSOA::PITCH*iy] = weno_plus(a[ix+L*iy], b[ix+L*iy], c[ix+L*iy], d[ix+L*iy], e[ix+L*iy],
                    b00, b01, b02, b10, b11, b12, b20, b21, b22,
                    c00, c01, c02, c10, c11, c12, c20, c21, c22,
                    d0, d1, d2);
}


/* *
 * Right hand side computation for the transport equations of the advected
 * quantities.  The computation of the velocity at the cell interface follows
 * E. Johnsen and T. Colonius, "Implementation of WENO schemes in compressible
 * multicomponent flow problems.", Journal of Computational Physics, (2006).
 * */
void Convection_CPP_HLLC_5eq_nonuniform::_xextraterm_hllc_nonuniform( const TempSOA& um, const TempSOA& up,
        const TempSOA& A2m, const TempSOA& A2p,
        const TempSOA& pm, const TempSOA& pp,
        const TempSOA& sm, const TempSOA& sp, const TempSOA& ss,
        const Real* const invh)
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

            divu.ref(ix, iy) = invh[ix]*(u_hllc1 - u_hllc0);
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

void Convection_CPP_HLLC_5eq_nonuniform::_yextraterm_hllc_nonuniform( const TempSOA& vm, const TempSOA& vp,
        const TempSOA& A2m, const TempSOA& A2p,
        const TempSOA& pm, const TempSOA& pp,
        const TempSOA& sm, const TempSOA& sp, const TempSOA& ss,
        const Real* const invh)
{
    for (int iy = 0; iy < OutputSOA::NY; ++iy)
    {
        const Real ihy = invh[iy];
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

            divu.ref(ix, iy) += ihy*(v_hllc1 - v_hllc0);
        }
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

void Convection_CPP_HLLC_5eq_nonuniform::_zextraterm_hllc_nonuniform( const TempSOA& wm0, const TempSOA& wp0,
        const TempSOA& wm1, const TempSOA& wp1,
        const TempSOA& A2m,  const TempSOA& A2p,
        const TempSOA& pm, const TempSOA& pp,
        const TempSOA& sm0, const TempSOA& sp0, const TempSOA& ss0,
        const TempSOA& sm1, const TempSOA& sp1, const TempSOA& ss1,
        const Real invh)
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

            divu.ref(ix, iy) += invh*(w_hllc1 - w_hllc0);
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


void Convection_CPP_HLLC_5eq_nonuniform::_xdivergence_nonuniform(const TempSOA& flux, OutputSOA& rhs, const Real* const invh)
{
    const Real * const f = flux.ptr(0,0);
    Real * const r = &rhs.ref(0,0);

    for(int iy=0; iy<OutputSOA::NY; iy++)
        for(int ix=0; ix<OutputSOA::NX; ix++)
            r[ix + OutputSOA::PITCH*iy] = invh[ix]*(f[ix + 1 + TempSOA::PITCH*iy] - f[ix + TempSOA::PITCH*iy]);
}

void Convection_CPP_HLLC_5eq_nonuniform::_ydivergence_nonuniform(const TempSOA& flux, OutputSOA& rhs, const Real* const invh)
{
    const Real * const f = flux.ptr(0,0);
    Real * const r = &rhs.ref(0,0);

    for(int iy=0; iy<OutputSOA::NY; iy++)
    {
        const Real ihy = invh[iy];
        for(int ix=0; ix<OutputSOA::NX; ix++)
            r[ix + OutputSOA::PITCH*iy] += ihy*(f[iy +  1 + TempSOA::PITCH*ix] - f[iy + TempSOA::PITCH*ix]);
    }
}

void Convection_CPP_HLLC_5eq_nonuniform::_zdivergence_nonuniform(const TempSOA& fback, const TempSOA& fforward, OutputSOA& rhs, const Real invh)
{
    const Real * const ff = fforward.ptr(0,0);
    const Real * const fb = fback.ptr(0,0);

    Real * const r = &rhs.ref(0,0);

    for(int iy=0; iy<OutputSOA::NY; iy++)
        for(int ix=0; ix<OutputSOA::NX; ix++)
            r[ix + OutputSOA::PITCH*iy] += invh*(ff[ix +  TempSOA::PITCH*iy] - fb[ix + TempSOA::PITCH*iy]);
}


void Convection_CPP_HLLC_5eq_nonuniform::_xflux_nonuniform(const int relid, const Coefficients_t& coeffs_minus, const Coefficients_t& coeffs_plus, const Real* const invh)
{
#ifndef _A2_A2R2_R_
    _xweno_minus(rho1.ring(relid), rho1.weno.ref(0), coeffs_minus);
    _xweno_pluss(rho1.ring(relid), rho1.weno.ref(1), coeffs_plus);
#endif
    _xweno_minus(rho2.ring(relid), rho2.weno.ref(0), coeffs_minus);
    _xweno_pluss(rho2.ring(relid), rho2.weno.ref(1), coeffs_plus);
#ifdef _A2_A2R2_R_
    _xweno_minus(rho.ring(relid), rho.weno.ref(0), coeffs_minus);
    _xweno_pluss(rho.ring(relid), rho.weno.ref(1), coeffs_plus);
#endif
    _xweno_minus(u.ring(relid), u.weno.ref(0), coeffs_minus);
    _xweno_pluss(u.ring(relid), u.weno.ref(1), coeffs_plus);
    _xweno_minus(v.ring(relid), v.weno.ref(0), coeffs_minus);
    _xweno_pluss(v.ring(relid), v.weno.ref(1), coeffs_plus);
    _xweno_minus(w.ring(relid), w.weno.ref(0), coeffs_minus);
    _xweno_pluss(w.ring(relid), w.weno.ref(1), coeffs_plus);
    _xweno_minus(p.ring(relid), p.weno.ref(0), coeffs_minus);
    _xweno_pluss(p.ring(relid), p.weno.ref(1), coeffs_plus);
    _xweno_minus(A2.ring(relid), A2.weno.ref(0), coeffs_minus);
    _xweno_pluss(A2.ring(relid), A2.weno.ref(1), coeffs_plus);

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
    _xextraterm_hllc_nonuniform(u.weno(0), u.weno(1), A2.weno(0), A2.weno(1), p.weno(0), p.weno(1), charvel(0), charvel(1), charvel_star(0), invh);

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


void Convection_CPP_HLLC_5eq_nonuniform::_yflux_nonuniform(const int relid, const Coefficients_t& coeffs_minus, const Coefficients_t& coeffs_plus, const Real* const invh)
{
#ifndef _A2_A2R2_R_
    _yweno_minus(rho1.ring(relid), rho1.weno.ref(0), coeffs_minus);
    _yweno_pluss(rho1.ring(relid), rho1.weno.ref(1), coeffs_plus);
#endif
    _yweno_minus(rho2.ring(relid), rho2.weno.ref(0), coeffs_minus);
    _yweno_pluss(rho2.ring(relid), rho2.weno.ref(1), coeffs_plus);
#ifdef _A2_A2R2_R_
    _yweno_minus(rho.ring(relid), rho.weno.ref(0), coeffs_minus);
    _yweno_pluss(rho.ring(relid), rho.weno.ref(1), coeffs_plus);
#endif
    _yweno_minus(u.ring(relid), u.weno.ref(0), coeffs_minus);
    _yweno_pluss(u.ring(relid), u.weno.ref(1), coeffs_plus);
    _yweno_minus(v.ring(relid), v.weno.ref(0), coeffs_minus);
    _yweno_pluss(v.ring(relid), v.weno.ref(1), coeffs_plus);
    _yweno_minus(w.ring(relid), w.weno.ref(0), coeffs_minus);
    _yweno_pluss(w.ring(relid), w.weno.ref(1), coeffs_plus);
    _yweno_minus(p.ring(relid), p.weno.ref(0), coeffs_minus);
    _yweno_pluss(p.ring(relid), p.weno.ref(1), coeffs_plus);
    _yweno_minus(A2.ring(relid), A2.weno.ref(0), coeffs_minus);
    _yweno_pluss(A2.ring(relid), A2.weno.ref(1), coeffs_plus);

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
    _yextraterm_hllc_nonuniform(v.weno(0), v.weno(1), A2.weno(0), A2.weno(1), p.weno(0), p.weno(1), charvel(0), charvel(1), charvel_star(0), invh);

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


void Convection_CPP_HLLC_5eq_nonuniform::_zflux_nonuniform(const int relid, const bool evalXtra, const int sliceid, const Coefficients_t& coeffs_minus, const Coefficients_t& coeffs_plus, const Real invh)
{
#ifndef _A2_A2R2_R_
    _zweno_minus(relid, rho1.ring, rho1.weno.ref(0), sliceid, coeffs_minus);
    _zweno_pluss(relid, rho1.ring, rho1.weno.ref(1), sliceid, coeffs_plus);
#endif
    _zweno_minus(relid, rho2.ring, rho2.weno.ref(0), sliceid, coeffs_minus);
    _zweno_pluss(relid, rho2.ring, rho2.weno.ref(1), sliceid, coeffs_plus);
#ifdef _A2_A2R2_R_
    _zweno_minus(relid, rho.ring, rho.weno.ref(0), sliceid, coeffs_minus);
    _zweno_pluss(relid, rho.ring, rho.weno.ref(1), sliceid, coeffs_plus);
#endif
    _zweno_minus(relid, u.ring, u.weno.ref(0), sliceid, coeffs_minus);
    _zweno_pluss(relid, u.ring, u.weno.ref(1), sliceid, coeffs_plus);
    _zweno_minus(relid, v.ring, v.weno.ref(0), sliceid, coeffs_minus);
    _zweno_pluss(relid, v.ring, v.weno.ref(1), sliceid, coeffs_plus);
    _zweno_minus(relid, w.ring, w.weno.ref(0), sliceid, coeffs_minus);
    _zweno_pluss(relid, w.ring, w.weno.ref(1), sliceid, coeffs_plus);
    _zweno_minus(relid, p.ring, p.weno.ref(0), sliceid, coeffs_minus);
    _zweno_pluss(relid, p.ring, p.weno.ref(1), sliceid, coeffs_plus);
    _zweno_minus(relid, A2.ring, A2.weno.ref(0), sliceid, coeffs_minus);
    _zweno_pluss(relid, A2.ring, A2.weno.ref(1), sliceid, coeffs_plus);

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
    if (evalXtra) _zextraterm_hllc_nonuniform(w.weno(-2), w.weno(-1), w.weno(0), w.weno(1), A2.weno(-1), A2.weno(0), p.weno(-1), p.weno(0), charvel(-2), charvel(-1), charvel_star(-1), charvel(0), charvel(1), charvel_star(0), invh);

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


void Convection_CPP_HLLC_5eq_nonuniform::_xrhs_nonuniform(const Real* const invh)
{
    _xdivergence_nonuniform(rho1.flux(), rho1.rhs, invh);
    _xdivergence_nonuniform(rho2.flux(), rho2.rhs, invh);
    _xdivergence_nonuniform(u.flux(), u.rhs, invh);
    _xdivergence_nonuniform(v.flux(), v.rhs, invh);
    _xdivergence_nonuniform(w.flux(), w.rhs, invh);
    _xdivergence_nonuniform(p.flux(), p.rhs, invh);
    _xdivergence_nonuniform(A2.flux(), A2.rhs, invh);
}

void Convection_CPP_HLLC_5eq_nonuniform::_yrhs_nonuniform(const Real* const invh)
{
    _ydivergence_nonuniform(rho1.flux(), rho1.rhs, invh);
    _ydivergence_nonuniform(rho2.flux(), rho2.rhs, invh);
    _ydivergence_nonuniform(u.flux(), u.rhs, invh);
    _ydivergence_nonuniform(v.flux(), v.rhs, invh);
    _ydivergence_nonuniform(w.flux(), w.rhs, invh);
    _ydivergence_nonuniform(p.flux(), p.rhs, invh);
    _ydivergence_nonuniform(A2.flux(), A2.rhs, invh);
}

void Convection_CPP_HLLC_5eq_nonuniform::_zrhs_nonuniform(const Real invh)
{
    _zdivergence_nonuniform(rho1.flux(-1), rho1.flux(0), rho1.rhs, invh);
    _zdivergence_nonuniform(rho2.flux(-1), rho2.flux(0), rho2.rhs, invh);
    _zdivergence_nonuniform(u.flux(-1), u.flux(0), u.rhs, invh);
    _zdivergence_nonuniform(v.flux(-1), v.flux(0), v.rhs, invh);
    _zdivergence_nonuniform(w.flux(-1), w.flux(0), w.rhs, invh);
    _zdivergence_nonuniform(p.flux(-1), p.flux(0), p.rhs, invh);
    _zdivergence_nonuniform(A2.flux(-1), A2.flux(0), A2.rhs, invh);
}


void Convection_CPP_HLLC_5eq_nonuniform::_copyback(Real * const gptfirst, const int gptfloats, const int rowgpts)
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

            rhs.r1 = a*rhs.r1 - m_dt*rho1.rhs(ix, iy);
            rhs.r2 = a*rhs.r2 - m_dt*rho2.rhs(ix, iy);
            rhs.u  = a*rhs.u -  m_dt*u.rhs(ix, iy);
            rhs.v  = a*rhs.v -  m_dt*v.rhs(ix, iy);
            rhs.w  = a*rhs.w -  m_dt*w.rhs(ix, iy);
            rhs.E  = a*rhs.E -  m_dt*p.rhs(ix, iy);
#ifdef _NOK_
            rhs.A2 = a*rhs.A2 - m_dt*(A2.rhs(ix, iy) - divu(ix,iy)*sumA2(ix,iy)*factor2);
#else
            rhs.A2 = a*rhs.A2 - m_dt*(A2.rhs(ix, iy) - divu(ix,iy)*sumA2(ix,iy)*factor2 - divu(ix,iy)*sumK(ix,iy)*factor2);
#endif
            // printf("%e\t%e\t%e\t%e\t%e\t%e\t%e\n", rhs.r1, rhs.r2, rhs.u, rhs.v, rhs.w, rhs.E, rhs.A2);
        }
}
