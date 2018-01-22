// Here we compute the smoothness indicators according to Jiang & Shu 1996.
// The difference between WenoFused_QPX.h and this implementation is that
// WenoFused_QPX.h factors the squared parenthesis into individual terms.
#pragma once

#include "common.h"

#ifdef _QPX_
#define vec_min(a, b)   vec_sel(a, b, vec_cmpgt(a, b))
#define vec_max(a, b)   vec_sel(a, b, vec_cmplt(a, b))
#endif

struct Weno_QPX_fused
{
    void weno_minus_plus_fused_opt2(const vector4double a, const vector4double b, const vector4double c, const vector4double d, const vector4double e, const vector4double f, Real *minus, Real *plus) //106 + 6*division + 2*reciprocal FLOP
    {
        // smoothness indicators
        ///////////////////////////////////////////////////////////////////////
        const vector4double F_1_4  = vec_splats(0.25f);
        const vector4double F_3    = vec_splats(3.0f);
        const vector4double F_13_3 = vec_splats(13.0f/3.0f);
        const vector4double M_2    = vec_splats(-2.0f);
        const vector4double M_4    = vec_splats(-4.0f);

        const vector4double mpt12_22 = vec_madd(c, M_2, vec_add(b, d));
        const vector4double mpt22_12 = vec_madd(d, M_2, vec_add(c, e));

        // minus
        const vector4double mt01 = vec_madd(c, F_3, vec_madd(b, M_4, a));
        const vector4double mt02 = vec_madd(b, M_2, vec_add(a, c));
        const vector4double m_is0 = vec_madd(F_13_3, vec_mul(mt02,mt02), vec_mul(mt01,mt01));

        const vector4double mt11 = vec_sub(b, d);
        const vector4double m_is1 = vec_madd(F_13_3, vec_mul(mpt12_22,mpt12_22), vec_mul(mt11,mt11));

        const vector4double mt21 = vec_madd(d, M_4, vec_madd(c, F_3, e));
        const vector4double m_is2 = vec_madd(F_13_3, vec_mul(mpt22_12,mpt22_12), vec_mul(mt21,mt21));

        // plus
        const vector4double pt01 = vec_madd(d, F_3, vec_madd(e, M_4, f));
        const vector4double pt02 = vec_madd(e, M_2, vec_add(f, d));
        const vector4double p_is0 = vec_madd(F_13_3, vec_mul(pt02,pt02), vec_mul(pt01,pt01));

        const vector4double pt11 = vec_sub(e, c);
        const vector4double p_is1 = vec_madd(F_13_3, vec_mul(mpt22_12,mpt22_12), vec_mul(pt11,pt11));

        const vector4double pt21 = vec_madd(c, M_4, vec_madd(d, F_3, b));
        const vector4double p_is2 = vec_madd(F_13_3, vec_mul(mpt12_22,mpt12_22), vec_mul(pt21,pt21));

        const vector4double WENOEPS4 = vec_splats(WENOEPS);

        // indicator
        const vector4double m_is0plus = vec_madd(F_1_4, m_is0, WENOEPS4);
        const vector4double m_is1plus = vec_madd(F_1_4, m_is1, WENOEPS4);
        const vector4double m_is2plus = vec_madd(F_1_4, m_is2, WENOEPS4);

        const vector4double p_is0plus = vec_madd(F_1_4, p_is0, WENOEPS4);
        const vector4double p_is1plus = vec_madd(F_1_4, p_is1, WENOEPS4);
        const vector4double p_is2plus = vec_madd(F_1_4, p_is2, WENOEPS4);
        ///////////////////////////////////////////////////////////////////////
        // 66

        // unnormalized non-linear weights
        ///////////////////////////////////////////////////////////////////////
        const vector4double m_alpha0 = mydivision<preclevel>(vec_splats(0.1f), vec_mul(m_is0plus, m_is0plus));
        const vector4double m_alpha1 = mydivision<preclevel>(vec_splats(0.6f), vec_mul(m_is1plus, m_is1plus));
        const vector4double m_alpha2 = mydivision<preclevel>(vec_splats(0.3f), vec_mul(m_is2plus, m_is2plus));
        const vector4double m_inv_alpha = myreciprocal<preclevel>(vec_add(m_alpha0, vec_add(m_alpha1, m_alpha2)));

        const vector4double p_alpha0 = mydivision<preclevel>(vec_splats(0.1f), vec_mul(p_is0plus, p_is0plus));
        const vector4double p_alpha1 = mydivision<preclevel>(vec_splats(0.6f), vec_mul(p_is1plus, p_is1plus));
        const vector4double p_alpha2 = mydivision<preclevel>(vec_splats(0.3f), vec_mul(p_is2plus, p_is2plus));
        const vector4double p_inv_alpha = myreciprocal<preclevel>(vec_add(p_alpha0, vec_add(p_alpha1, p_alpha2)));
        ///////////////////////////////////////////////////////////////////////
        // 10 + 6*division + 2*reciprocal

        {
            const vector4double F_2_6  = vec_splats(2.0f/6.0f);
            const vector4double F_5_2  = vec_splats(2.5f);
            const vector4double F_11_2 = vec_splats(5.5f);
            const vector4double M_1_2  = vec_splats(-0.5f);
            const vector4double M_7_2  = vec_splats(-3.5f);

            // clipping
            ///////////////////////////////////////////////////////////////////////
            const vector4double min_val_minus = vec_min(a, vec_min(b, vec_min(c, vec_min(d,e))));
            const vector4double max_val_minus = vec_max(a, vec_max(b, vec_max(c, vec_max(d,e))));
            const vector4double min_val_plus = vec_min(f, vec_min(b, vec_min(c, vec_min(d,e))));
            const vector4double max_val_plus = vec_max(f, vec_max(b, vec_max(c, vec_max(d,e))));
            ///////////////////////////////////////////////////////////////////////

            const vector4double m1_p2 = vec_madd(b, M_1_2, vec_madd(c, F_5_2, d));
            const vector4double m2_p1 = vec_madd(e, M_1_2, vec_madd(d, F_5_2, c));
            const vector4double m0 = vec_madd(b, M_7_2, vec_madd(c, F_11_2, a));
            const vector4double p0 = vec_madd(e, M_7_2, vec_madd(d, F_11_2, f));

            const vector4double minus_tmp = vec_madd(m_alpha0, m0, vec_madd(m_alpha1, m1_p2, vec_mul(m_alpha2, m2_p1)));
            const vector4double plus_tmp  = vec_madd(p_alpha0, p0, vec_madd(p_alpha1, m2_p1, vec_mul(p_alpha2, m1_p2)));

            vec_sta(vec_max(min_val_minus, vec_min(max_val_minus, vec_mul(minus_tmp, vec_mul(F_2_6, m_inv_alpha)))), 0L, minus);
            vec_sta(vec_max(min_val_plus, vec_min(max_val_plus, vec_mul(plus_tmp , vec_mul(F_2_6, p_inv_alpha)))), 0L, plus);
            // 30
        }
    }
};
