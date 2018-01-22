/*
 *  Source_CPP.h
 *  MPCFcore
 *
 *  Created by Ursula Rasthofer
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#pragma once
#include <cassert>
#include <cstdlib>
#include <cstdio>

#include "common.h"
#include "BlockInfo.h"

class Source_CPP
{
protected:

    Real a, dt;
    Real g[3], f[3];

    inline bool _is_aligned(const void * const ptr, unsigned int alignment) const
    {
        return ((size_t)ptr) % alignment == 0;
    }

public:

    Source_CPP(Real a, Real dt, const Real gaccel[3], const Real volforce[3]): a(a), dt(dt)
    {
        g[0]=gaccel[0];
        g[1]=gaccel[1];
        g[2]=gaccel[2];
        f[0]=volforce[0];
        f[1]=volforce[1];
        f[2]=volforce[2];
    }

    void compute (const Real * const srcfirst, const int srcfloats, Real * const dstfirst, const int dstfloats) const;

    static void printflops(const float PEAKPERF_CORE, const float PEAKBAND, const size_t NCORES, const size_t NT, const size_t NBLOCKS, const float MEASUREDTIME, const bool bAwk=false)
    {}
};


class OneWayAcousticSource_CPP
{
public:

    struct SourceParameter
    {
        Real gamma, smooth;
        Real t0, c0, x0, x0a;
        Real amplitude, amplitudea, omega, T;
        Real sigma_t;
        Real phase0, phase0a, frequency, ambient;
    };

    OneWayAcousticSource_CPP(const Real t, const Real dt, const SourceParameter& p):
        t(t), dt(dt), param(p)
    {}

    void compute(const Real * const srcfirst, const int srcfloats, Real * const dstfirst, const int dstfloats, const BlockInfo& info) const
    {
        assert(srcfloats==dstfloats);

        const Real f = _f(t);
        const Real fa = _fa(t);

        const Real reg0 = 1. / (std::sqrt(2. * M_PI) * param.smooth);
        const Real reg1 = 1. / (2. * std::pow(param.smooth, 2));

        for (int iz = 0; iz < _BLOCKSIZE_; ++iz)
            for (int iy = 0; iy < _BLOCKSIZE_; ++iy)
                for (int ix = 0; ix < _BLOCKSIZE_; ++ix)
                {
                    const int i = dstfloats * (ix + _BLOCKSIZE_ * (iy + _BLOCKSIZE_ * iz));
                    Real pos[3];
                    info.pos(pos, ix, iy, iz);

                    //std::cout << "*** " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
                    //std::terminate();

                    // if (std::abs(pos[0] - param.x0) < 6.0*param.smooth)
                    {
                        const Real reg = reg0 * std::exp(-std::pow(pos[0] - param.x0, 2) * reg1);
                        const Real rega = reg0 * std::exp(-std::pow(pos[0] - param.x0a, 2) * reg1);
                        /*
                        if (iz == 0 && iy == 0 && ix == 0) {
                          std::cout << "***" 
                            << " x=" << pos[0] 
                            << " y=" << pos[1] 
                            << " z=" << pos[2] 
                            << " reg=" << reg 
                            << " x0=" << param.x0 
                            << " smooth=" << param.smooth 
                            << " reg0=" << reg0
                            << " reg1=" << reg1 
                            << " f=" << f 
                            << " freq=" << param.frequency
                            << " t=" << t 
                            << " t0=" << param.t0
                            << " sigma_t=" << param.sigma_t
                            << std::endl;
                        }
                        */
                        dstfirst[i+0] += dt * f * reg / param.c0; // mass 1
                        dstfirst[i+2] += dt * f * reg; // x-momentum
                        dstfirst[i+5] += dt * f * reg * param.c0 / (param.gamma - 1.); // energy

                        dstfirst[i+0] += dt * fa * rega / param.c0; // mass 1
                        dstfirst[i+2] += -dt * fa * rega; // x-momentum
                        dstfirst[i+5] += dt * fa * rega * param.c0 / (param.gamma - 1.); // energy
                    }
                }
    }

protected:

    const Real t;
    const Real dt;
    SourceParameter param;

    inline bool _is_aligned(const void * const ptr, unsigned int alignment) const
    {
        return ((size_t)ptr) % alignment == 0;
    }

    // source function
    inline Real _f(const Real t) const
    {
        return param.amplitude 
          * std::exp(-std::pow((t - param.t0) * param.frequency / param.sigma_t, 2))
          * std::cos(2.0 * M_PI * (param.frequency * (t - param.t0) - param.phase0));
    }
    inline Real _fa(const Real t) const
    {
        return param.amplitudea 
          * std::exp(-std::pow((t - param.t0) * param.frequency / param.sigma_t, 2))
          * std::cos(2.0 * M_PI * (param.frequency * (t - param.t0) - param.phase0a));
    }
};
