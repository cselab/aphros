/*
 *  Diffusion_CPP.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 2/27/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */
#pragma once
#include <cmath>

#include "DivTensor_CPP.h"

class Diffusion_CPP: public virtual DivTensor_CPP
{

// currently unused as unnecessary (rasthofer August 2016)
/*
    void _convert_single_phase(const Real * const gptfirst, const int gptfloats, const int rowgpts);
    void _convert_two_phase(const Real * const gptfirst, const int gptfloats, const int rowgpts);
*/

protected:

    struct AssumedType { Real r1, r2, u, v, w, E, A2, dummy; };

    typedef Real RealTemp;

// see comment above
/*
    inline void _convert(const Real * const gptfirst, const int gptfloats, const int rowgpts)
    {
        if (mu1 == mu2) _convert_single_phase(gptfirst, gptfloats, rowgpts);
        else            _convert_two_phase(gptfirst, gptfloats, rowgpts);
    }
*/
    void _convert(const Real * const gptfirst, const int gptfloats, const int rowgpts);

    virtual void _xface(const TempSOA_ST& ux0, const TempSOA_ST& uy0, const TempSOA_ST& uz0, const TempSOA_ST& ux1, const TempSOA_ST& uy1, const TempSOA_ST& uz1,
            const TempSOA_ST& vx0, const TempSOA_ST& vy0, const TempSOA_ST& vx1, const TempSOA_ST& vy1,
            const TempSOA_ST& wx0, const TempSOA_ST& wz0, const TempSOA_ST& wx1, const TempSOA_ST& wz1,
            const InputSOA_ST &mu0);

    virtual void _yface(const TempSOA_ST& ux0, const TempSOA_ST& uy0, const TempSOA_ST& ux1, const TempSOA_ST& uy1,
            const TempSOA_ST& vx0, const TempSOA_ST& vy0, const TempSOA_ST& vz0, const TempSOA_ST& vx1, const TempSOA_ST& vy1, const TempSOA_ST& vz1,
            const TempSOA_ST& wy0, const TempSOA_ST& wz0, const TempSOA_ST& wy1, const TempSOA_ST& wz1,
            const InputSOA_ST &mu0);

    virtual void _zface(const TempSOA_ST& ux, const TempSOA_ST& uz,
            const TempSOA_ST& vy, const TempSOA_ST& vz,
            const TempSOA_ST& wx, const TempSOA_ST& wy, const TempSOA_ST& wz,
            const InputSOA_ST &mu0, const InputSOA_ST &mu1,
            TempPiZSOA_ST& tzx, TempPiZSOA_ST& tzy, TempPiZSOA_ST& tzz);

    inline void _grad_next()
    {
        ringnx.next(); ringny.next(); ringnz.next();
        gradvx.next(); gradvy.next(); gradvz.next();
        gradwx.next(); gradwy.next(); gradwz.next();
    }

public:

    const RealTemp mu1, mu2;
    const Real smoothing_length;

    RingTempSOA_ST gradvx, gradvy, gradvz, gradwx, gradwy, gradwz;

    Diffusion_CPP(const Real a,
            const Real mu1,
            const Real mu2,
            const Real h,
            const Real smoothing_length,
            const Real dtinvh):
        DivTensor_CPP(a, dtinvh, h, 1.0), mu1(mu1), mu2(mu2), smoothing_length(smoothing_length) { }

    void compute(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
            Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts)
    {
        _convert(srcfirst, srcfloats, rowsrcs);
        _input_next();

        _convert(srcfirst + srcfloats*slicesrcs, srcfloats, rowsrcs);

        _corners(ringu(-1), ringu(), ringnx.ref(), ringny.ref(), ringnz.ref());
        _corners(ringv(-1), ringv(), gradvx.ref(), gradvy.ref(), gradvz.ref());
        _corners(ringw(-1), ringw(), gradwx.ref(), gradwy.ref(), gradwz.ref());

        _zface(ringnx(), ringnz(), gradvy(), gradvz(),
                gradwx(), gradwy(), gradwz(),
                ringls(-1), ringls(),
                ringtzx.ref(), ringtzy.ref(), ringtzz.ref());

        _udot_tz(ringu(-1), ringv(-1), ringw(-1),  ringu(), ringv(), ringw(), ringtzx(), ringtzy(), ringtzz(), ringutz.ref());

        for(int islice=0; islice<_BLOCKSIZE_; islice++)
        {
            _tensors_next();
            _input_next();
            _grad_next();

            _convert(srcfirst + (islice+2)*srcfloats*slicesrcs, srcfloats, rowsrcs);

            _corners(ringu(-1), ringu(0), ringnx.ref(), ringny.ref(), ringnz.ref());
            _corners(ringv(-1), ringv(0), gradvx.ref(), gradvy.ref(), gradvz.ref());
            _corners(ringw(-1), ringw(0), gradwx.ref(), gradwy.ref(), gradwz.ref());

            _xface(ringnx(-1), ringny(-1), ringnz(-1), ringnx(), ringny(), ringnz(),
                    gradvx(-1), gradvy(-1), gradvx(), gradvy(),
                    gradwx(-1), gradwz(-1), gradwx(), gradwz(), ringls(-1));
            /* _xface(ringnx(-1), ringny(-1), ringnz(-1), ringnx(), ringny(), ringnz(), */
            /*         gradvx(-1), gradvy(-1), gradvx(), gradvy(), */
            /*         gradwx(-1), gradwz(-1), gradwx(), gradwz()); */

            _yface(ringnx(-1), ringny(-1), ringnx(), ringny(),
                    gradvx(-1), gradvy(-1), gradvz(-1), gradvx(), gradvy(), gradvz(),
                    gradwy(-1), gradwz(-1), gradwy(), gradwz(), ringls(-1));
            /* _yface(ringnx(-1), ringny(-1), ringnx(), ringny(), */
            /*         gradvx(-1), gradvy(-1), gradvz(-1), gradvx(), gradvy(), gradvz(), */
            /*         gradwy(-1), gradwz(-1), gradwy(), gradwz()); */

            _zface(ringnx(), ringnz(), gradvy(), gradvz(),
                    gradwx(), gradwy(), gradwz(),
                    ringls(-1), ringls(),
                    ringtzx.ref(), ringtzy.ref(), ringtzz.ref());
            /* _zface(ringnx(), ringnz(), gradvy(), gradvz(), */
            /*         gradwx(), gradwy(), gradwz(), ringtzx.ref(), ringtzy.ref(), ringtzz.ref()); */

            _udot_tx(ringu(-1), ringv(-1), ringw(-1));
            _udot_ty(ringu(-1), ringv(-1), ringw(-1));
            _udot_tz(ringu(-1), ringv(-1), ringw(-1), ringu(0), ringv(0), ringw(0), ringtzx(), ringtzy(), ringtzz(), ringutz.ref());

            _div_dxy();
            _div_dz(ringtzx(-1), ringtzy(-1), ringtzz(-1), ringutz(-1), ringtzx(0), ringtzy(0), ringtzz(0), ringutz(0));

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
        DivTensor_CPP::hpc_info(flop_convert, traffic_convert, flop_corners, traffic_corner, flop_tensorface, traffic_tensorface,
                flop_tdotu, traffic_tdotu, flop_div, traffic_div, flop_copyback, traffic_copyback, footprint);

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
