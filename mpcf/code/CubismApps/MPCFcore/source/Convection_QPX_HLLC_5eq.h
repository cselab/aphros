/*
 *  Convection_QPX_HLLC_5eq.h
 *  MPCFcore
 *
 *  Created by Fabian Wermelinger 03/02/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef CONVECTION_QPX_HLLC_5EQ_H_3IIUORPM
#define CONVECTION_QPX_HLLC_5EQ_H_3IIUORPM

#include "check_errors.h"
#include "Convection_CPP_HLLC_5eq.h"
#include "WenoSOA2D_QPX.h"
#include "HLLCSOA2D_QPX_5eq.h"
#include "DivSOA2D_QPX.h"

class Convection_QPX_HLLC_5eq : public Convection_CPP_HLLC_5eq
{
protected:

    #if defined(_BGQ_)
        __align(_ALIGNBYTES_) struct TinyScratchPad { Real tmp[4][4];};
        #else
        struct __align(_ALIGNBYTES_) TinyScratchPad { Real tmp[4][4];};
    #endif

    void _qpx_convert_aligned(Real * const gptfirst, const int gptfloats, const int rowgpts,
                              Real * const ar1, Real * const ar2, Real * const u, Real * const v,
                              Real * const w, Real * const p, Real * const a2)
    {
        const vector4double F_1 = vec_splats(1);
        const vector4double F_0 = vec_splats(0);

#define DESTID (dx + (InputSOA::PITCH)*dy)

        for(int dy=0; dy<_BLOCKSIZE_+6; dy++)
        {
            Real * const in = gptfirst + dy*gptfloats*rowgpts -gptfloats;

            for(int dx=0; dx<_BLOCKSIZE_+8; dx+=4)
            {
                const int WID = (dx + (int)(dx==0))*gptfloats; // west
                const int CID = (dx+1)*gptfloats; // center
                const int EID = (dx+3 - (int)(dx==_BLOCKSIZE_+4))*gptfloats; // east

                vector4double dataA0 = vec_lda(0L, in + WID);
                vector4double dataA1 = vec_lda(0L, in + WID + 4);
                vector4double dataB0 = vec_lda(0L, in + CID);
                vector4double dataB1 = vec_lda(0L, in + CID + 4);
                vector4double dataC0 = vec_lda(0L, in + CID + gptfloats);
                vector4double dataC1 = vec_lda(0L, in + CID + gptfloats + 4);
                vector4double dataD0 = vec_lda(0L, in + EID);
                vector4double dataD1 = vec_lda(0L, in + EID + 4);

                _DIEGO_TRANSPOSE4(dataA0, dataB0, dataC0, dataD0);

#ifdef _CONVERTCLIP_
                dataA0 = mymax(F_0, dataA0);
                dataB0 = mymax(F_0, dataB0);
#endif /* _CONVERTCLIP_ */

                vec_sta(dataA0, 0L, ar1 + DESTID);
                vec_sta(dataB0, 0L, ar2 + DESTID);

                const vector4double inv_rho = myreciprocal<preclevel>(vec_add(dataA0,dataB0));

                vec_sta(vec_mul(dataC0, inv_rho), 0L, u + DESTID);
                vec_sta(vec_mul(dataD0, inv_rho), 0L, v + DESTID);

                _DIEGO_TRANSPOSE4(dataA1, dataB1, dataC1, dataD1);

#if defined(_ALPHACLIP_) || defined(_CONVERTCLIP_)
                // dataC1 = mymax(F_0, mymin(F_1, dataC1)); // alpha_2
                dataC1 = mymax(vec_splats((Real)ALPHAEPS), mymin(vec_splats(1.0f-(Real)ALPHAEPS), dataC1)); // alpha2
#endif /* _CONVERTCLIP_ */

                const vector4double alpha1 = vec_sub(F_1, dataC1); // alpha_1

                const vector4double GmixInv = myreciprocal<preclevel>(vec_madd(alpha1, vec_splats(_g1m1Inv), vec_mul(dataC1, vec_splats(_g2m1Inv))));

#ifdef _CONVERTCLIP_
                const vector4double Pmix = mymax(F_0, vec_madd(alpha1, vec_splats(_g1m1Inv*_g1*_pc1), vec_mul(dataC1, vec_splats(_g2m1Inv*_g2*_pc2))));
#else
                const vector4double Pmix = vec_madd(alpha1, vec_splats(_g1m1Inv*_g1*_pc1), vec_mul(dataC1, vec_splats(_g2m1Inv*_g2*_pc2)));
#endif /* _CONVERTCLIP_ */

                const vector4double fac = vec_mul(vec_splats(-0.5f), inv_rho);

                vec_sta(vec_mul(dataA1, inv_rho), 0L, w + DESTID);

                const vector4double speedsquared = vec_madd(dataC0, dataC0, vec_madd(dataD0, dataD0, vec_mul(dataA1, dataA1)));
                const vector4double actpres = vec_mul(GmixInv, vec_sub(vec_madd(fac, speedsquared, dataB1), Pmix));

#ifdef _CONVERTCLIP_
                const vector4double fab = vec_sub(vec_splats(PRESEPS), actpres);
#ifdef _NOK_
                const vector4double pthresh = vec_mul(vec_mul(Pmix, GmixInv), myreciprocal<preclevel>(vec_add(GmixInv, F_1)));
#else
                const vector4double stage1 = vec_sel(mymin(vec_splats(_pc1), vec_splats(_pc2)), vec_splats(_pc2), vec_sub(dataC1, F_1));
                const vector4double pthresh = vec_sel(vec_splats(_pc1), stage1, vec_cmpgt(dataC1, vec_splats(0.0f)));
#endif
                const vector4double pdiff = vec_sub(fab, pthresh);
                const vector4double pclip = vec_add(actpres, vec_sel(vec_splats(0.0f), pdiff, pdiff));
#endif /* _CONVERTCLIP_ */

#ifdef _CONVERTCLIP_
                vec_sta(pclip, 0L, p + DESTID);
#else
                vec_sta(actpres, 0L, p + DESTID);
#endif /* _CONVERTCLIP_ */
                vec_sta(dataC1, 0L, a2 + DESTID);
            }
        }

#undef DESTID

    }

    virtual void _convert(const Real * const gptfirst, const int gptfloats, const int rowgpts, const int islice)
    {
        {
            const size_t x = (size_t)gptfirst;
            const int remainder = x & 0x0f;
            assert(remainder == 0);
            assert(gptfloats == 8);
            if (remainder)
            {
                printf("oooops! pointer is not aligned: 0x%x\n", (int)x);
                abort();
            }
        }

        _qpx_convert_aligned(const_cast<Real*>(gptfirst), gptfloats, rowgpts,
                             & rho1.ring.ref().ref(-4,-3),
                             & rho2.ring.ref().ref(-4,-3),
                             & u.ring.ref().ref(-4,-3),
                             & v.ring.ref().ref(-4,-3),
                             & w.ring.ref().ref(-4,-3),
                             & p.ring.ref().ref(-4,-3),
                             & A2.ring.ref().ref(-4,-3));

    }

    virtual void _xflux(const int relid)
    {
        {
            WenoSOA2D_QPX wenoizer;

            wenoizer.xcompute(rho1.ring(relid), rho1.weno.ref(0), rho1.weno.ref(1));
            wenoizer.xcompute(rho2.ring(relid), rho2.weno.ref(0), rho2.weno.ref(1));
            wenoizer.xcompute(u.ring(relid), u.weno.ref(0), u.weno.ref(1));
            wenoizer.xcompute(v.ring(relid), v.weno.ref(0), v.weno.ref(1));
            wenoizer.xcompute(w.ring(relid), w.weno.ref(0), w.weno.ref(1));
            wenoizer.xcompute(p.ring(relid), p.weno.ref(0), p.weno.ref(1));
            wenoizer.xcompute(A2.ring(relid), A2.weno.ref(0), A2.weno.ref(1));
        }

        HLLCSOA2D_QPX_5eq hllczator(_g1,_g2,_pc1,_pc2);
        hllczator.all(rho1.weno(0), rho1.weno(1),
                rho2.weno(0), rho2.weno(1),
                u.weno(0), u.weno(1),
                v.weno(0), v.weno(1),
                w.weno(0), w.weno(1),
                p.weno(0), p.weno(1),
                A2.weno(0), A2.weno(1),
                charvel.ref(0), charvel.ref(1), charvel_star.ref(0),
                rho1.flux.ref(), rho2.flux.ref(),
                u.flux.ref(), v.flux.ref(), w.flux.ref(), p.flux.ref(),
                A2.flux.ref());

        DivSOA2D_HLLC_QPX_5eq divtor(_g1,_g2,_pc1,_pc2);
        divtor.xextraterm(u.weno(0), u.weno(1), A2.weno(0), A2.weno(1), p.weno(0), p.weno(1), charvel(0), charvel(1), charvel_star(0), divu, sumA2, sumK);
    }

    virtual void _yflux(const int relid)
    {
        {
            WenoSOA2D_QPX wenoizer;

            wenoizer.ycompute(rho1.ring(relid), rho1.weno.ref(0), rho1.weno.ref(1));
            wenoizer.ycompute(rho2.ring(relid), rho2.weno.ref(0), rho2.weno.ref(1));
            wenoizer.ycompute(u.ring(relid), u.weno.ref(0), u.weno.ref(1));
            wenoizer.ycompute(v.ring(relid), v.weno.ref(0), v.weno.ref(1));
            wenoizer.ycompute(w.ring(relid), w.weno.ref(0), w.weno.ref(1));
            wenoizer.ycompute(p.ring(relid), p.weno.ref(0), p.weno.ref(1));
            wenoizer.ycompute(A2.ring(relid), A2.weno.ref(0), A2.weno.ref(1));
        }

        HLLCSOA2D_QPX_5eq hllczator(_g1,_g2,_pc1,_pc2);
        hllczator.all(rho1.weno(0), rho1.weno(1),
                rho2.weno(0), rho2.weno(1),
                v.weno(0), v.weno(1),
                u.weno(0), u.weno(1),
                w.weno(0), w.weno(1),
                p.weno(0), p.weno(1),
                A2.weno(0), A2.weno(1),
                charvel.ref(0), charvel.ref(1), charvel_star.ref(0),
                rho1.flux.ref(), rho2.flux.ref(),
                v.flux.ref(), u.flux.ref(), w.flux.ref(), p.flux.ref(),
                A2.flux.ref());

        DivSOA2D_HLLC_QPX_5eq divtor(_g1,_g2,_pc1,_pc2);
        divtor.yextraterm(v.weno(0), v.weno(1), A2.weno(0), A2.weno(1), p.weno(0), p.weno(1), charvel(0), charvel(1), charvel_star(0), divu, sumA2, sumK);
    }

    virtual void _zflux(const int relid, const bool evalXtra)
    {
        {
            WenoSOA2D_QPX wenoizer;

            wenoizer.zcompute(relid, rho1.ring, rho1.weno.ref(0), rho1.weno.ref(1));
            wenoizer.zcompute(relid, rho2.ring, rho2.weno.ref(0), rho2.weno.ref(1));
            wenoizer.zcompute(relid, u.ring, u.weno.ref(0), u.weno.ref(1));
            wenoizer.zcompute(relid, v.ring, v.weno.ref(0), v.weno.ref(1));
            wenoizer.zcompute(relid, w.ring, w.weno.ref(0), w.weno.ref(1));
            wenoizer.zcompute(relid, p.ring, p.weno.ref(0), p.weno.ref(1));
            wenoizer.zcompute(relid, A2.ring, A2.weno.ref(0), A2.weno.ref(1));
        }

        HLLCSOA2D_QPX_5eq hllczator(_g1,_g2,_pc1,_pc2);
        hllczator.all(rho1.weno(0), rho1.weno(1),
                rho2.weno(0), rho2.weno(1),
                w.weno(0), w.weno(1),
                u.weno(0), u.weno(1),
                v.weno(0), v.weno(1),
                p.weno(0), p.weno(1),
                A2.weno(0), A2.weno(1),
                charvel.ref(0), charvel.ref(1), charvel_star.ref(0),
                rho1.flux.ref(), rho2.flux.ref(),
                w.flux.ref(), u.flux.ref(), v.flux.ref(), p.flux.ref(),
                A2.flux.ref());

        DivSOA2D_HLLC_QPX_5eq divtor(_g1,_g2,_pc1,_pc2);
        if (evalXtra)
            divtor.zextraterm(w.weno(-2), w.weno(-1),
                    w.weno(0), w.weno(1),
                    A2.weno(-1), A2.weno(0),
                    p.weno(-1), p.weno(0),
                    charvel(-2), charvel(-1),
                    charvel(0), charvel(1),
                    charvel_star(-1), charvel_star(0),
                    divu, sumA2, sumK);
    }

    virtual void _xrhs()
    {
        DivSOA2D_QPX divtor;
        divtor.xrhs(rho1.flux(), rho1.rhs);
        divtor.xrhs(rho2.flux(), rho2.rhs);
        divtor.xrhs(u.flux(), u.rhs);
        divtor.xrhs(v.flux(), v.rhs);
        divtor.xrhs(w.flux(), w.rhs);
        divtor.xrhs(p.flux(), p.rhs);
        divtor.xrhs(A2.flux(), A2.rhs);
    }

    virtual void _yrhs()
    {
        DivSOA2D_QPX divtor;
        divtor.yrhs(rho1.flux(), rho1.rhs);
        divtor.yrhs(rho2.flux(), rho2.rhs);
        divtor.yrhs(u.flux(), u.rhs);
        divtor.yrhs(v.flux(), v.rhs);
        divtor.yrhs(w.flux(), w.rhs);
        divtor.yrhs(p.flux(), p.rhs);
        divtor.yrhs(A2.flux(), A2.rhs);
    }

    virtual void _zrhs()
    {
        DivSOA2D_QPX divtor;
        divtor.zrhs(rho1.flux(-1), rho1.flux(0), rho1.rhs);
        divtor.zrhs(rho2.flux(-1), rho2.flux(0), rho2.rhs);
        divtor.zrhs(u.flux(-1), u.flux(0), u.rhs);
        divtor.zrhs(v.flux(-1), v.flux(0), v.rhs);
        divtor.zrhs(w.flux(-1), w.flux(0), w.rhs);
        divtor.zrhs(p.flux(-1), p.flux(0), p.rhs);
        divtor.zrhs(A2.flux(-1), A2.flux(0), A2.rhs);
    }

    virtual void _copyback(Real * const gptfirst, const int gptfloats, const int rowgpts)
    {
        const vector4double mya = vec_splats(a);
        const vector4double lambda = vec_splats(dtinvh);
        const vector4double M_1_6 = vec_splats(-1.f/6);

        const int offset1 = gptfloats;
        const int offset2 = 2 * gptfloats;
        const int offset3 = 3 * gptfloats;

        for(int iy=0; iy<OutputSOA::NY; iy++)
        {
            Real * const r1ptr = &rho1.rhs.ref(0, iy);
            Real * const r2ptr = &rho2.rhs.ref(0, iy);
            Real * const uptr = &u.rhs.ref(0, iy);
            Real * const vptr = &v.rhs.ref(0, iy);
            Real * const wptr = &w.rhs.ref(0, iy);
            Real * const pptr = &p.rhs.ref(0, iy);
            Real * const A2ptr = &A2.rhs.ref(0, iy);
            Real * const sumA2ptr = &sumA2.ref(0, iy);
#ifndef _NOK_
            Real * const sumKptr = &sumK.ref(0, iy);
#endif /* _NOK_ */
            Real * const divuptr = &divu.ref(0, iy);

            for(int ix=0; ix<OutputSOA::NX; ix += 4)
            {
                const int entry_out = gptfloats*(ix + iy*rowgpts);

                vector4double d0 = vec_mul(lambda, vec_lda(0L, r1ptr + ix));
                vector4double d1 = vec_mul(lambda, vec_lda(0L, r2ptr + ix));
                vector4double d2 = vec_mul(lambda, vec_lda(0L, uptr + ix));
                vector4double d3 = vec_mul(lambda, vec_lda(0L, vptr + ix));
                vector4double d0B = vec_mul(lambda, vec_lda(0L, wptr + ix));
                vector4double d1B = vec_mul(lambda, vec_lda(0L, pptr + ix));
                vector4double d2B = vec_mul(M_1_6, vec_lda(0L, sumA2ptr + ix));
                vector4double d3B = vec_splats(0);

                vector4double mydivu = vec_lda(0L, divuptr + ix);
#ifdef _NOK_
                d2B = vec_mul(lambda, vec_madd(d2B, mydivu, vec_lda(0L, A2ptr + ix)));
#else
                const vector4double sK = vec_mul(M_1_6, vec_lda(0L, sumKptr + ix));
                d2B = vec_mul(lambda, vec_madd(sK, mydivu, vec_madd(d2B, mydivu, vec_lda(0L, A2ptr + ix))));
#endif /* _NOK_ */

                _DIEGO_TRANSPOSE4(d0, d1, d2, d3);
                _DIEGO_TRANSPOSE4(d0B, d1B, d2B, d3B);

                d0  = vec_msub(mya, vec_lda(0L, gptfirst + entry_out), d0);
                d0B = vec_msub(mya, vec_lda(4 * sizeof(Real), gptfirst + entry_out), d0B);
                d1  = vec_msub(mya, vec_lda(0L, gptfirst + entry_out + offset1), d1);
                d1B = vec_msub(mya, vec_lda(4 * sizeof(Real), gptfirst + entry_out + offset1), d1B);
                d2  = vec_msub(mya, vec_lda(0L, gptfirst + entry_out + offset2), d2);
                d2B = vec_msub(mya, vec_lda(4 * sizeof(Real), gptfirst + entry_out + offset2), d2B);
                d3  = vec_msub(mya, vec_lda(0L, gptfirst + entry_out + offset3), d3);
                d3B = vec_msub(mya, vec_lda(4 * sizeof(Real), gptfirst + entry_out + offset3), d3B);

                vec_sta(d0, 0L, gptfirst + entry_out);
                vec_sta(d0B, 4 * sizeof(Real), gptfirst + entry_out);
                vec_sta(d1, 0L, gptfirst + entry_out + offset1);
                vec_sta(d1B, 4 * sizeof(Real), gptfirst + entry_out + offset1);
                vec_sta(d2, 0L, gptfirst + entry_out + offset2);
                vec_sta(d2B, 4 * sizeof(Real), gptfirst + entry_out + offset2);
                vec_sta(d3, 0L, gptfirst + entry_out + offset3);
                vec_sta(d3B, 4 * sizeof(Real), gptfirst + entry_out + offset3);
            }
        }
    }

public:
    Convection_QPX_HLLC_5eq(const Real a, const Real dtinvh, const Real g1=1.4, const Real g2=1.4, const Real pc1=0.0, const Real pc2=0.0) :
    Convection_CPP_HLLC_5eq(a, dtinvh, g1, g2, pc1, pc2) {}
};

#endif /* CONVECTION_QPX_HLLC_5EQ_H_3IIUORPM */
