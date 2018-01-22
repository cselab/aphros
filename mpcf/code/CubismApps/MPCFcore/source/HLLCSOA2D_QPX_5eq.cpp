/*
 *  HLLCSOA2D_QPX_5eq.cpp
 *  MPCFcore
 *
 *  Created by Fabian Wermelinger 03/03/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */

#include "HLLCSOA2D_QPX_5eq.h"

using namespace std;

#ifndef _MICROFUSION_
#define _MICROFUSION_ 2
#endif

inline void _qpx_hllc_all_5eq(Real * const r1m, Real * const r1p,
        Real * const r2m, Real * const r2p,
        Real * const vdm, Real * const vdp,
        Real * const v1m, Real * const v1p,
        Real * const v2m, Real * const v2p,
        Real * const pm, Real * const pp,
        Real * const A2m, Real * const A2p,
        Real * const outam, Real * const outap, Real * const outas,
        Real * const outrho1, Real * const outrho2,
        Real * const outvd, Real * const outv1, Real * const outv2,
        Real * const oute, Real * const outA2,
        const Real g1, const Real pc1, const Real g2, const Real pc2)
{
    enum { NTOTAL = TempSOA::PITCH * TempSOA::NY};

    const vector4double F_0 = vec_splats(0.0f);
    const vector4double F_1 = vec_splats(1.0f);

    /* const vector4double G1 = vec_splats(g1); */
    /* const vector4double G2 = vec_splats(g2); */
    /* const vector4double P1 = vec_splats(pc1); */
    /* const vector4double P2 = vec_splats(pc2); */

    const vector4double GI1 = vec_splats(1.0f/(g1-1.0f));
    const vector4double GI2 = vec_splats(1.0f/(g2-1.0f));


    for(int ID = 0; ID < NTOTAL; ID += 4)
    {
        // ********************************************************************
        // MINUS PREP
        // ********************************************************************
#ifdef _RECONPCLIP_ // implies _ALPHACLIP_
        // const vector4double A2minus = mymax(F_0, mymin(F_1, vec_lda(0L, A2m + ID)));
        const vector4double A2minus = mymax(vec_splats((Real)ALPHAEPS), mymin(vec_splats(1.0f-(Real)ALPHAEPS), vec_lda(0L, A2m + ID)));
        vector4double pminusTMP = vec_lda(0L, pm + ID);
#else
#ifdef _ALPHACLIP_
        const vector4double A2minus = mymax(vec_splats((Real)ALPHAEPS), mymin(vec_splats(1.0f-(Real)ALPHAEPS), vec_lda(0L, A2m + ID)));
#else
        const vector4double A2minus = vec_lda(0L, A2m + ID);
#endif /* _ALPHACLIP_ */
        const vector4double pminus = vec_lda(0L, pm + ID);
#endif /* _RECONPCLIP_ */

        const vector4double Gmixminus = vec_madd(vec_sub(F_1, A2minus), GI1, vec_mul(A2minus, GI2));

#ifdef _RECONPCLIP_
        const vector4double Pmixminus = mymax(F_0, vec_madd(vec_sub(F_1,A2minus), vec_mul(vec_mul(GI1,vec_splats(g1)),vec_splats(pc1)), vec_mul(A2minus,vec_mul(vec_mul(GI2,vec_splats(g2)),vec_splats(pc2)))));
        const vector4double fabminus = vec_sub(vec_splats(PRESEPS), pminusTMP);
    #ifdef _NOK_
        const vector4double GmixInvminus = myreciprocal<preclevel>(Gmixminus);
        const vector4double pthreshminus = vec_mul(vec_mul(Pmixminus, GmixInvminus), myreciprocal<preclevel>(vec_add(GmixInvminus, F_1)));
    #else
        const vector4double stage1minus = vec_sel(mymin(vec_splats(pc1), vec_splats(pc2)), vec_splats(pc2), vec_sub(A2minus, F_1));
        const vector4double pthreshminus = vec_sel(vec_splats(pc1), stage1minus, vec_cmpgt(A2minus, F_0));
    #endif /* _NOK_ */
        const vector4double pdiffminus = vec_sub(fabminus, pthreshminus);
        const vector4double pminus = vec_add(pminusTMP, vec_sel(F_0, pdiffminus, pdiffminus));
#else
    #ifdef _NOK_
        const vector4double GmixInvminus = myreciprocal<preclevel>(Gmixminus);
    #endif /* _NOK_ */
        const vector4double Pmixminus = vec_madd(vec_sub(F_1,A2minus), vec_mul(vec_mul(GI1,vec_splats(g1)),vec_splats(pc1)), vec_mul(A2minus,vec_mul(vec_mul(GI2,vec_splats(g2)),vec_splats(pc2))));
#endif /* _RECONPCLIP_ */

        // load minus
        const vector4double r1minus = vec_lda(0L, r1m + ID);
        const vector4double r2minus = vec_lda(0L, r2m + ID);
        const vector4double vdminus = vec_lda(0L, vdm + ID);
        const vector4double v1minus = vec_lda(0L, v1m + ID);
        const vector4double v2minus = vec_lda(0L, v2m + ID);

        // some intermediates minus
        const vector4double rminus     = vec_add(r1minus, r2minus);
        const vector4double speedminus = vec_madd(vdminus, vdminus, vec_madd(v1minus, v1minus, vec_mul(v2minus, v2minus)));
        const vector4double uminus     = vec_mul(vdminus, rminus);
        const vector4double uminus_v1  = vec_mul(v1minus, rminus);
        const vector4double uminus_v2  = vec_mul(v2minus, rminus);
        const vector4double eminus     = vec_madd(pminus, Gmixminus, vec_madd(vec_mul(vec_splats(0.5f), rminus), speedminus, Pmixminus));

        // Speed of sound
#ifdef _NOK_
        const vector4double cminus2    = vec_mul(vec_madd(GmixInvminus, Pmixminus, vec_madd(GmixInvminus, pminus, pminus)), myreciprocal<preclevel>(rminus));
#else

        const vector4double gpc1minus  = vec_mul(vec_splats(g1), vec_add(pminus, vec_splats(pc1)));
        const vector4double gpc2minus  = vec_mul(vec_splats(g2), vec_add(pminus, vec_splats(pc2)));
#if 0
        // VARIANT 1: (Similar as in MaxSpeedOfSound_QPX_5eq)
        const vector4double cminus2    = vec_mul(vec_mul(gpc1minus,gpc2minus), myreciprocal<preclevel>(vec_mul(rminus, vec_madd(vec_sub(F_1,A2minus), gpc2minus, vec_mul(A2minus,gpc1minus)))));
#else
        // VARIANT 2: (Similar as in MaxSpeedOfSound_QPX_5eq)
        const vector4double sel_1minus = vec_cmpgt(gpc1minus,F_0);
        const vector4double sel_2minus = vec_cmpgt(gpc2minus,F_0);

        const vector4double rc2minus_1 = vec_sel(F_1, gpc1minus, sel_1minus);
        const vector4double rc2minus_2 = vec_sel(F_1, gpc2minus, sel_2minus);
        const vector4double a_1minus   = vec_sel(F_0, vec_sub(F_1,A2minus), sel_1minus);
        const vector4double a_2minus   = vec_sel(F_0, A2minus, sel_2minus);

        const vector4double cminus2    = vec_mul(vec_mul(rc2minus_1,rc2minus_2), myreciprocal<preclevel>(vec_mul(rminus,vec_madd(a_1minus,rc2minus_2,vec_mul(a_2minus,rc2minus_1)))));
#endif

#endif /* _NOK_ */

        const vector4double cminus     = mysqrt<preclevel>(cminus2);
        // ********************************************************************


        // ********************************************************************
        // PLUS PREP
        // ********************************************************************
#ifdef _RECONPCLIP_ // implies _ALPHACLIP_
        // const vector4double A2plus = mymax(F_0, mymin(F_1, vec_lda(0L, A2p + ID)));
        const vector4double A2plus = mymax(vec_splats((Real)ALPHAEPS), mymin(vec_splats(1.0f-(Real)ALPHAEPS), vec_lda(0L, A2p + ID)));
        vector4double pplusTMP = vec_lda(0L, pp + ID);
#else
#ifdef _ALPHACLIP_
        const vector4double A2plus = mymax(vec_splats((Real)ALPHAEPS), mymin(vec_splats(1.0f-(Real)ALPHAEPS), vec_lda(0L, A2p + ID)));
#else
        const vector4double A2plus = vec_lda(0L, A2p + ID);
#endif /* _ALPHACLIP_ */
        const vector4double pplus = vec_lda(0L, pp + ID);
#endif /* _RECONPCLIP_ */

        const vector4double Gmixplus = vec_madd(vec_sub(F_1, A2plus), GI1, vec_mul(A2plus, GI2));

#ifdef _RECONPCLIP_
        const vector4double Pmixplus = mymax(F_0, vec_madd(vec_sub(F_1,A2plus), vec_mul(vec_mul(GI1,vec_splats(g1)),vec_splats(pc1)), vec_mul(A2plus,vec_mul(vec_mul(GI2,vec_splats(g2)),vec_splats(pc2)))));
        const vector4double fabplus = vec_sub(vec_splats(PRESEPS), pplusTMP);
    #ifdef _NOK_
        const vector4double GmixInvplus = myreciprocal<preclevel>(Gmixplus);
        const vector4double pthreshplus = vec_mul(vec_mul(Pmixplus, GmixInvplus), myreciprocal<preclevel>(vec_add(GmixInvplus, F_1)));
    #else
        const vector4double stage1plus = vec_sel(mymin(vec_splats(pc1), vec_splats(pc2)), vec_splats(pc2), vec_sub(A2plus, F_1));
        const vector4double pthreshplus = vec_sel(vec_splats(pc1), stage1plus, vec_cmpgt(A2plus, F_0));
    #endif /* _NOK_ */
        const vector4double pdiffplus = vec_sub(fabplus, pthreshplus);
        const vector4double pplus = vec_add(pplusTMP, vec_sel(F_0, pdiffplus, pdiffplus));
#else
    #ifdef _NOK_
        const vector4double GmixInvplus = myreciprocal<preclevel>(Gmixplus);
    #endif /* _NOK_ */
        const vector4double Pmixplus = vec_madd(vec_sub(F_1,A2plus), vec_mul(vec_mul(GI1,vec_splats(g1)),vec_splats(pc1)), vec_mul(A2plus,vec_mul(vec_mul(GI2,vec_splats(g2)),vec_splats(pc2))));
#endif /* _RECONPCLIP_ */

        // load plus
        const vector4double r1plus = vec_lda(0L, r1p + ID);
        const vector4double r2plus = vec_lda(0L, r2p + ID);
        const vector4double vdplus = vec_lda(0L, vdp + ID);
        const vector4double v1plus = vec_lda(0L, v1p + ID);
        const vector4double v2plus = vec_lda(0L, v2p + ID);

        // some intermediates plus
        const vector4double rplus     = vec_add(r1plus, r2plus);
        const vector4double speedplus = vec_madd(vdplus, vdplus, vec_madd(v1plus, v1plus, vec_mul(v2plus, v2plus)));
        const vector4double uplus     = vec_mul(vdplus, rplus);
        const vector4double uplus_v1  = vec_mul(v1plus, rplus);
        const vector4double uplus_v2  = vec_mul(v2plus, rplus);
        const vector4double eplus     = vec_madd(pplus, Gmixplus, vec_madd(vec_mul(vec_splats(0.5f), rplus), speedplus, Pmixplus));

        // Speed of sound
#ifdef _NOK_
        const vector4double cplus2    = vec_mul(vec_madd(GmixInvplus, Pmixplus, vec_madd(GmixInvplus, pplus, pplus)), myreciprocal<preclevel>(rplus));
#else

        const vector4double gpc1plus  = vec_mul(vec_splats(g1), vec_add(pplus, vec_splats(pc1)));
        const vector4double gpc2plus  = vec_mul(vec_splats(g2), vec_add(pplus, vec_splats(pc2)));
#if 0
        // VARIANT 1: (Similar as in MaxSpeedOfSound_QPX_5eq)
        const vector4double cplus2    = vec_mul(vec_mul(gpc1plus,gpc2plus), myreciprocal<preclevel>(vec_mul(rplus, vec_madd(vec_sub(F_1,A2plus), gpc2plus, vec_mul(A2plus,gpc1plus)))));
#else
        // VARIANT 2: (Similar as in MaxSpeedOfSound_QPX_5eq)
        const vector4double sel_1plus = vec_cmpgt(gpc1plus,F_0);
        const vector4double sel_2plus = vec_cmpgt(gpc2plus,F_0);

        const vector4double rc2plus_1 = vec_sel(F_1, gpc1plus, sel_1plus);
        const vector4double rc2plus_2 = vec_sel(F_1, gpc2plus, sel_2plus);
        const vector4double a_1plus   = vec_sel(F_0, vec_sub(F_1,A2plus), sel_1plus);
        const vector4double a_2plus   = vec_sel(F_0, A2plus, sel_2plus);

        const vector4double cplus2    = vec_mul(vec_mul(rc2plus_1,rc2plus_2), myreciprocal<preclevel>(vec_mul(rplus,vec_madd(a_1plus,rc2plus_2,vec_mul(a_2plus,rc2plus_1)))));
#endif

#endif /* _NOK_ */

        const vector4double cplus     = mysqrt<preclevel>(cplus2);
        // ********************************************************************


        // ********************************************************************
        // WAVESPEED ESTIMATES
        // ********************************************************************
        const vector4double Rr    = mysqrt<preclevel>(vec_mul(rplus, myreciprocal<preclevel>(rminus)));
        const vector4double Rinv  = myreciprocal<preclevel>(vec_add(F_1, Rr));
        const vector4double eta_2 = vec_mul(vec_mul(vec_mul(vec_splats(0.5f),Rr),Rinv),Rinv);
        const vector4double d2    = vec_madd(eta_2,vec_mul(vec_sub(vdplus,vdminus),vec_sub(vdplus, vdminus)), vec_mul(vec_madd(Rr,cplus2,cminus2),Rinv));
        const vector4double d     = mysqrt<preclevel>(d2);
        const vector4double u     = vec_mul(vec_madd(Rr,vdplus,vdminus),Rinv);

        const vector4double aminus = mymin(vec_sub(vdminus, cminus), vec_sub(u, d));
        const vector4double aplus = mymax(vec_add(u, d), vec_add(vdplus, cplus));
        vec_sta(aminus, 0L, outam + ID);
        vec_sta(aplus, 0L, outap + ID);

        const vector4double facm = vec_mul(rminus,vec_sub(aminus, vdminus));
        const vector4double facp = vec_mul(rplus,vec_sub(aplus, vdplus));
        const vector4double astar = vec_mul(vec_madd(vdminus,facm, vec_sub(vec_sub(pplus,pminus), vec_mul(vdplus,facp))), myreciprocal<preclevel>(vec_sub(facm,facp))) ;
        vec_sta(astar, 0L, outas + ID);


        // ********************************************************************
        // HLLC FLUXES
        // ********************************************************************
        /* *
         * The flux computation is split into 4 parts:
         * 1.) Compute signum of s^*, compute s^- and s^+
         * 2.) Compute chi^* and delta of q^* and q
         * 3.) Compute trivial flux
         * 4.) Compute HLLC flux
         * */

        // 1.)
        /* (fabianw; Thu 03 Mar 2016 07:45:33 PM CET) This one is a
         * cheaper version.  It neglects the case astar == 0 (might happen
         * rarely in general, but will compute a wrong sign if the case
         * happens) */
        /* TODO: (fabianw; Fri 04 Mar 2016 05:36:22 PM CET) is the return type
         * of vec_cmpgt float or integral? implicit cast here?? */
        // const vector4double sign_star = vec_cmpgt(astar, F_0);

        /* (fabianw; Thu 03 Mar 2016 07:45:51 PM CET) Correct version */
        const vector4double sign_star= vec_sel(vec_cmpgt(astar,F_0), F_0, vec_cmpeq(vec_abs(astar),F_0));

        const vector4double s_minus  = mymin(aminus,F_0);
        const vector4double s_pluss  = mymax(aplus, F_0);
        /* const vector4double s_minus  = vec_sel(aminus, F_0, aminus); */
        /* const vector4double s_pluss  = vec_sel(F_0, aplus, aplus); */


        // 2.)
        const vector4double chi_starm = vec_mul(vec_sub(aminus, vdminus), myreciprocal<preclevel>(vec_sub(aminus, astar)));
        const vector4double chi_starp = vec_mul(vec_sub(aplus, vdplus), myreciprocal<preclevel>(vec_sub(aplus, astar)));
        const vector4double q_deltam_rho1  = vec_msub(r1minus, chi_starm, r1minus);
        const vector4double q_deltap_rho1  = vec_msub(r1plus, chi_starp, r1plus);
        const vector4double q_deltam_rho2  = vec_msub(r2minus, chi_starm, r2minus);
        const vector4double q_deltap_rho2  = vec_msub(r2plus, chi_starp, r2plus);
        const vector4double q_deltam_vd  = vec_msub(rminus, vec_mul(chi_starm, astar), uminus);
        const vector4double q_deltap_vd  = vec_msub(rplus, vec_mul(chi_starp, astar), uplus);
        const vector4double q_deltam_v1  = vec_msub(uminus_v1, chi_starm, uminus_v1);
        const vector4double q_deltap_v1  = vec_msub(uplus_v1, chi_starp, uplus_v1);
        const vector4double q_deltam_v2  = vec_msub(uminus_v2, chi_starm, uminus_v2);
        const vector4double q_deltap_v2  = vec_msub(uplus_v2, chi_starp, uplus_v2);
        const vector4double q_deltam_e  = vec_msub(vec_madd(vec_madd(rminus, astar, vec_mul(pminus, myreciprocal<preclevel>(vec_sub(aminus, vdminus)))), vec_sub(astar,vdminus), eminus), chi_starm, eminus);
        const vector4double q_deltap_e  = vec_msub(vec_madd(vec_madd(rplus, astar, vec_mul(pplus, myreciprocal<preclevel>(vec_sub(aplus, vdplus)))), vec_sub(astar,vdplus), eplus), chi_starp, eplus);
        const vector4double q_deltam_A2  = vec_msub(A2minus, chi_starm, A2minus);
        const vector4double q_deltap_A2  = vec_msub(A2plus, chi_starp, A2plus);


        // 3.)
        const vector4double fminus_rho1 = vec_mul(vdminus, r1minus);
        const vector4double fminus_rho2 = vec_mul(vdminus, r2minus);
        const vector4double fminus_vd = vec_madd(vdminus, uminus, pminus);
        const vector4double fminus_v1 = vec_mul(vdminus, uminus_v1);
        const vector4double fminus_v2 = vec_mul(vdminus, uminus_v2);
        const vector4double fminus_e = vec_mul(vdminus, vec_add(pminus, eminus));
        const vector4double fminus_A2 = vec_mul(vdminus, A2minus);

        const vector4double fplus_rho1 = vec_mul(vdplus, r1plus);
        const vector4double fplus_rho2 = vec_mul(vdplus, r2plus);
        const vector4double fplus_vd = vec_madd(vdplus, uplus, pplus);
        const vector4double fplus_v1 = vec_mul(vdplus, uplus_v1);
        const vector4double fplus_v2 = vec_mul(vdplus, uplus_v2);
        const vector4double fplus_e = vec_mul(vdplus, vec_add(pplus, eplus));
        const vector4double fplus_A2 = vec_mul(vdplus, A2plus);


        // 4.)
        /* const vector4double factorm = vec_mul(vec_splats(0.5f),vec_add(F_1, sign_star)); */
        /* const vector4double factorp = vec_mul(vec_splats(0.5f),vec_sub(F_1, sign_star)); */
        const vector4double factorm = vec_madd(vec_splats(0.5f), sign_star, vec_splats(0.5f));
        const vector4double factorp = vec_madd(vec_splats(-0.5f), sign_star, vec_splats(0.5f));
        const vector4double result_rho1 = vec_add(vec_mul(factorm,vec_madd(s_minus,q_deltam_rho1,fminus_rho1)),vec_mul(factorp,vec_madd(s_pluss,q_deltap_rho1,fplus_rho1)));
        const vector4double result_rho2 = vec_add(vec_mul(factorm,vec_madd(s_minus,q_deltam_rho2,fminus_rho2)),vec_mul(factorp,vec_madd(s_pluss,q_deltap_rho2,fplus_rho2)));
        const vector4double result_vd = vec_add(vec_mul(factorm,vec_madd(s_minus,q_deltam_vd,fminus_vd)),vec_mul(factorp,vec_madd(s_pluss,q_deltap_vd,fplus_vd)));
        const vector4double result_v1 = vec_add(vec_mul(factorm,vec_madd(s_minus,q_deltam_v1,fminus_v1)),vec_mul(factorp,vec_madd(s_pluss,q_deltap_v1,fplus_v1)));
        const vector4double result_v2 = vec_add(vec_mul(factorm,vec_madd(s_minus,q_deltam_v2,fminus_v2)),vec_mul(factorp,vec_madd(s_pluss,q_deltap_v2,fplus_v2)));
        const vector4double result_e = vec_add(vec_mul(factorm,vec_madd(s_minus,q_deltam_e,fminus_e)),vec_mul(factorp,vec_madd(s_pluss,q_deltap_e,fplus_e)));
        const vector4double result_A2 = vec_add(vec_mul(factorm,vec_madd(s_minus,q_deltam_A2,fminus_A2)),vec_mul(factorp,vec_madd(s_pluss,q_deltap_A2,fplus_A2)));

        vec_sta(result_rho1, 0L, outrho1 + ID);
        vec_sta(result_rho2, 0L, outrho2 + ID);
        vec_sta(result_vd, 0L, outvd + ID);
        vec_sta(result_v1, 0L, outv1 + ID);
        vec_sta(result_v2, 0L, outv2 + ID);
        vec_sta(result_e, 0L, oute + ID);
        vec_sta(result_A2, 0L, outA2 + ID);
    }
}

void HLLCSOA2D_QPX_5eq::all(const TempSOA& r1minus, const TempSOA& r1plus,
                            const TempSOA& r2minus, const TempSOA& r2plus,
                            const TempSOA& vdminus, const TempSOA& vdplus,
                            const TempSOA& v1minus, const TempSOA& v1plus,
                            const TempSOA& v2minus, const TempSOA& v2plus,
                            const TempSOA& pminus, const TempSOA& pplus,
                            const TempSOA& A2minus, const TempSOA& A2plus,
                            TempSOA& outam, TempSOA& outap, TempSOA& outas,
                            TempSOA& outrho1, TempSOA& outrho2,
                            TempSOA& outvd, TempSOA& outv1, TempSOA& outv2,
                            TempSOA& oute,TempSOA& outA2) const
{
    _qpx_hllc_all_5eq(const_cast<Real*>(r1minus.ptr(0,0)), const_cast<Real*>(r1plus.ptr(0,0)),
            const_cast<Real*>(r2minus.ptr(0,0)), const_cast<Real*>(r2plus.ptr(0,0)),
            const_cast<Real*>(vdminus.ptr(0,0)), const_cast<Real*>(vdplus.ptr(0,0)),
            const_cast<Real*>(v1minus.ptr(0,0)), const_cast<Real*>(v1plus.ptr(0,0)),
            const_cast<Real*>(v2minus.ptr(0,0)), const_cast<Real*>(v2plus.ptr(0,0)),
            const_cast<Real*>(pminus.ptr(0,0)), const_cast<Real*>(pplus.ptr(0,0)),
            const_cast<Real*>(A2minus.ptr(0,0)), const_cast<Real*>(A2plus.ptr(0,0)),
            &outam.ref(0,0), &outap.ref(0,0), &outas.ref(0,0),
            &outrho1.ref(0,0), &outrho2.ref(0,0),
            &outvd.ref(0,0), &outv1.ref(0,0),  &outv2.ref(0,0),
            &oute.ref(0,0), &outA2.ref(0,0),
            m_g1, m_pc1, m_g2, m_pc2);
}
