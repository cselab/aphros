/*
 *  MaxSpeedOfSound_QPX_5eq.h
 *  MPCFcore
 *
 *  Created by Fabian Wermelinger 03/21/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef MAXSPEEDOFSOUND_QPX_5EQ_H_A8MWCLIV
#define MAXSPEEDOFSOUND_QPX_5EQ_H_A8MWCLIV

#include "MaxSpeedOfSound_CPP_5eq.h"

class MaxSpeedOfSound_QPX_5eq : public MaxSpeedOfSound_CPP_5eq
{
    enum
    {
        NPOINTS = _BLOCKSIZE_ * _BLOCKSIZE_ * _BLOCKSIZE_,
        JUMP = sizeof(Real) * 4
    };

    inline vector4double _compute_p(const vector4double invr, const vector4double E,
                                    const vector4double GmixInv, const vector4double Pmix,
                                    const vector4double speed2, const vector4double alpha2) const
    {
        const vector4double fac = vec_mul(vec_splats(-0.5f), invr);
        const vector4double actpres = vec_mul(GmixInv, vec_sub(vec_madd(fac, speed2, E), Pmix));

#ifdef _CONVERTCLIP_
        const vector4double fab = vec_sub(vec_splats(PRESEPS), actpres);
#ifdef _NOK_
        const vector4double pthresh = vec_mul(vec_mul(Pmix, GmixInv), myreciprocal<preclevel>(vec_add(GmixInv, vec_splats(1.0f))));
#else
        // TODO: (fabianw@mavt.ethz.ch; Mon 02 May 2016 02:44:05 PM CEST) what
        // to do here if ALPHAEPS is used?
        const vector4double stage1 = vec_sel(mymin(vec_splats(_pc1), vec_splats(_pc2)), vec_splats(_pc2), vec_sub(alpha2, vec_splats(1.0f)));
        const vector4double pthresh = vec_sel(vec_splats(_pc1), stage1, vec_cmpgt(alpha2, vec_splats(0.0f)));
#endif
        const vector4double pdiff = vec_sub(fab, pthresh);
#endif /* _CONVERTCLIP_ */

#ifdef _CONVERTCLIP_
        return vec_add(actpres, vec_sel(vec_splats(0.0f), pdiff, pdiff));
#else
        return actpres;
#endif /* _CONVERTCLIP_ */

    }

    inline vector4double _compute_c(const vector4double invr, const vector4double p,
                                    const vector4double GmixInv, const vector4double Pmix,
                                    const vector4double alpha2) const
    {
// #if 1
#ifdef _NOK_
        const vector4double c2 = vec_mul(vec_madd(GmixInv, Pmix, vec_madd(GmixInv, p, p)), invr);
#else
        // TODO: (fabianw@mavt.ethz.ch; Wed 17 Aug 2016 06:51:50 PM CEST) The
        // maximum speed of sound is used to determine the time step based on
        // the CFL number.  We choose a safety margin in the CFL number anyway,
        // hence for the computation here preclevel=0 may be sufficient in any
        // case
#if 0
        // VARIANT 1: (roughly 1.3X slower than _NOK_, roughly 1.2X faster than VARIANT 2)
        // The following 3 lines just divide using gpc1 and gpc2. No checks for
        // division by zero are made:
        const vector4double gpc1= vec_mul(vec_splats(_g1), vec_add(p, vec_splats(_pc1)));
        const vector4double gpc2= vec_mul(vec_splats(_g2), vec_add(p, vec_splats(_pc2)));
        const vector4double c2 = vec_mul(vec_mul(gpc1,gpc2), vec_mul(invr, myreciprocal<preclevel>(vec_madd(vec_sub(vec_splats(1.0f),alpha2), gpc2, vec_mul(alpha2,gpc1)))));
#else
        // VARIANT 2: (roughly 1.5X slower than _NOK_)
        // Checks wether gpc1 or gpc2 is zero.  This prevents NaN's but note
        // that if gpc1 or gpc2 is zero or less than zero, something else went
        // wrong already.
        const vector4double F_0 = vec_splats(0.0f);
        const vector4double F_1 = vec_splats(1.0f);

        const vector4double gpc1= vec_mul(vec_splats(_g1), vec_add(p, vec_splats(_pc1)));
        const vector4double gpc2= vec_mul(vec_splats(_g2), vec_add(p, vec_splats(_pc2)));

        const vector4double sel_1 = vec_cmpgt(gpc1,F_0);
        const vector4double sel_2 = vec_cmpgt(gpc2,F_0);
        const vector4double rc2_1 = vec_sel(F_1, gpc1, sel_1);
        const vector4double rc2_2 = vec_sel(F_1, gpc2, sel_2);
        const vector4double a_1   = vec_sel(F_0, vec_sub(F_1,alpha2), sel_1);
        const vector4double a_2   = vec_sel(F_0, alpha2, sel_2);

        const vector4double c2 = vec_mul(vec_mul(rc2_1,rc2_2), vec_mul(invr, myreciprocal<preclevel>(vec_madd(a_1,rc2_2,vec_mul(a_2,rc2_1)))));
#endif
#endif /* _NOK_ */

        return mysqrt<preclevel>(c2);
    }

    template<int GPTFLOATS> inline vector4double _sweep4(Real * const src) const
    {
        enum { PTJUMP = sizeof(Real) * GPTFLOATS };

        vector4double data0 = vec_lda(0L, src);         //alpha1rho1
        vector4double data1 = vec_lda(PTJUMP, src);     //alpha2rho2
        vector4double data2 = vec_lda(PTJUMP * 2, src); //ru
        vector4double data3 = vec_lda(PTJUMP * 3, src); //rv

        _DIEGO_TRANSPOSE4(data0, data1, data2, data3);

#ifdef _CONVERTCLIP_
        data0 = mymax(vec_splats(0.0f), data0);
        data1 = mymax(vec_splats(0.0f), data1);
#endif /* _CONVERTCLIP_ */

        const vector4double invr = myreciprocal<preclevel>(vec_add(data0,data1));
        vector4double speed2 = vec_madd(data2, data2, vec_mul(data3, data3));
        vector4double maxvel = mymax(vec_abs(data2), vec_abs(data3));

        data0 = vec_lda(JUMP, src);                     //rw
        data1 = vec_lda(PTJUMP + JUMP, src);            //E
        data2 = vec_lda(PTJUMP * 2 + JUMP, src);        //alpha2
        data3 = vec_lda(PTJUMP * 3 + JUMP, src);        //mickey mouse

        _DIEGO_TRANSPOSE4(data0, data1, data2, data3);

#if defined(_ALPHACLIP_) || defined(_CONVERTCLIP_)
        // data2 = mymax(vec_splats(0.0f), mymin(vec_splats(1.0f), data2)); // alpha_2
        data2 = mymax(vec_splats((Real)ALPHAEPS), mymin(vec_splats(1.0f-(Real)ALPHAEPS), data2)); // alpha_2
#endif /* _CONVERTCLIP_ */

        speed2 = vec_madd(data0, data0, speed2);
        maxvel = mymax(vec_abs(data0), maxvel);

        const vector4double alpha1 = vec_sub(vec_splats(1.0f), data2); // alpha_1
        const vector4double GmixInv = myreciprocal<preclevel>(vec_madd(alpha1, vec_splats(_g1m1Inv), vec_mul(data2, vec_splats(_g2m1Inv))));

#ifdef _CONVERTCLIP_
        const vector4double Pmix = mymax(vec_splats(0.0f), vec_madd(alpha1, vec_splats(_g1m1Inv*_g1*_pc1), vec_mul(data2, vec_splats(_g2m1Inv*_g2*_pc2))));
#else
        const vector4double Pmix = vec_madd(alpha1, vec_splats(_g1m1Inv*_g1*_pc1), vec_mul(data2, vec_splats(_g2m1Inv*_g2*_pc2)));
#endif /* _CONVERTCLIP_ */

        const vector4double p = _compute_p(invr, data1, GmixInv, Pmix, speed2, data2);
        const vector4double c = _compute_c(invr, p, GmixInv, Pmix, data2);

        return vec_madd(maxvel, invr, c);
    }

    template<int GPTFLOATS> Real _compute(Real * const src) const
    {
        enum { NFLOATS = GPTFLOATS * NPOINTS };

        vector4double sos4 = vec_splats(0);

        for(int i=0; i < NFLOATS; i += 4 * GPTFLOATS)
            sos4 = mymax(sos4, _sweep4<GPTFLOATS>(src + i));

        sos4 = mymax(sos4, vec_perm(sos4, sos4, vec_gpci(2323)));

        /* why does this not work?
         * sos4 = mymax(sos4, vec_perm(sos4, sos4, vec_gpci(1111)));
         * return vec_extract(sos4, 0);
         */

        return std::max(vec_extract(sos4, 0), vec_extract(sos4, 1));
    }


public:
    MaxSpeedOfSound_QPX_5eq(const Real g1=1.4, const Real g2=1.4, const Real pc1=0.0, const Real pc2=0.0) :
        MaxSpeedOfSound_CPP_5eq(g1,g2,pc1,pc2) {}

    Real compute(const Real * const _src, const int gptfloats) const
    {
        assert(gptfloats == 8 || gptfloats == 16);

        Real * const src = const_cast<Real *>(_src);

        if (gptfloats == 8)
            return _compute<8>(src);
        else if (gptfloats == 16)
            return _compute<16>(src);
        else
        {
            printf("ooops MaxSpeedOfSound_QPX_5eq::compute: gptfloats is not quite right. aborting.\n");
            abort();
            return -1;
        }
    }
};

#endif /* MAXSPEEDOFSOUND_QPX_5EQ_H_A8MWCLIV */
