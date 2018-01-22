/*
 *  Weno_QPX_3rdOrder.h
 *
 *
 *  Created by Diego Rossinelli on 4/12/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#ifdef _QPXEMU_
#include "QPXEMU.h"
#endif

#ifdef _QPX_
#define vec_min(a, b)   vec_sel(a, b, vec_cmpgt(a, b))
#define vec_max(a, b)   vec_sel(a, b, vec_cmplt(a, b))
#endif

struct WenoQPX_3rdOrder_Minus
{
protected:

	vector4double mywenoeps, F_1_3, F_2_3, F_1_2, F_3;

public:
	WenoQPX_3rdOrder_Minus()
	{
		mywenoeps = vec_splats(WENOEPS);

		F_1_3 = vec_splats(1.f/3);
		F_2_3 = vec_splats(2.f/3);
		F_1_2 = vec_splats(1.f/2);
		F_3 = vec_splats(3.f);
	}

	inline vector4double _weno3(const vector4double b, const vector4double c, const vector4double d) const
	{
        vector4double min_val = vec_min(b, vec_min(c, d));
        vector4double max_val = vec_max(b, vec_max(c, d));

		const vector4double cmb = vec_sub(c, b);
		const vector4double is0_eps = vec_madd(cmb, cmb, mywenoeps);

		const vector4double dmc = vec_sub(d, c);
		const vector4double is1_eps = vec_madd(dmc, dmc, mywenoeps);

		const vector4double a0 = vec_mul(F_1_3, myreciprocal<preclevel>(vec_mul(is0_eps, is0_eps)));
		const vector4double a1 = vec_mul(F_2_3, myreciprocal<preclevel>(vec_mul(is1_eps, is1_eps)));
		const vector4double inv_asum = myreciprocal<preclevel>(vec_add(a0, a1));

		const vector4double term0 = vec_msub(F_3, c, b);
		const vector4double term1 = vec_add(c, d);

        return vec_max(min_val, vec_min(max_val, vec_mul(vec_mul(F_1_2, inv_asum), vec_madd(a0, term0, vec_mul(a1, term1)))));
	}

	inline vector4double _weno(const vector4double a, const vector4double b, const vector4double c,
							   const vector4double d, const vector4double e) const
	{
		return _weno3(b, c, d);
	}
};

struct WenoQPX_3rdOrder_Plus : WenoQPX_3rdOrder_Minus
{
	WenoQPX_3rdOrder_Plus(): WenoQPX_3rdOrder_Minus() { }

	inline vector4double _weno(const vector4double a, const vector4double b, const vector4double c,
							   const vector4double d, const vector4double e) const
	{
		return _weno3(e, d, c);
	}
};

struct Weno_QPX_fused : WenoQPX_3rdOrder_Minus
{
	Weno_QPX_fused() : WenoQPX_3rdOrder_Minus() { }

	void weno_minus_plus_fused_opt2(const vector4double a, const vector4double b,
									const vector4double c, const vector4double d,
									const vector4double e, const vector4double f,
									Real *minus, Real *plus) const
	{
        vector4double min_val_minus = vec_min(b, vec_min(c, d));
        vector4double max_val_minus = vec_max(b, vec_max(c, d));

        vector4double min_val_plus = vec_min(c, vec_min(d, e));
        vector4double max_val_plus = vec_max(c, vec_max(d, e));

		//naive would be:
		//vec_sta(_weno3(b, c, d), 0L, minus);
		//vec_sta(_weno3(e, d, c), 0L, plus);

		const vector4double dmc = vec_sub(d, c);
		const vector4double is1_eps = vec_madd(dmc, dmc, mywenoeps);
		const vector4double a1 = vec_mul(F_2_3, myreciprocal<preclevel>(vec_mul(is1_eps, is1_eps)));
		const vector4double term1 = vec_add(c, d);

		//minus. nothing to share with plus here.
		{
			const vector4double cmb = vec_sub(c, b);
			const vector4double is0_eps = vec_madd(cmb, cmb, mywenoeps);
			const vector4double a0 = vec_mul(F_1_3, myreciprocal<preclevel>(vec_mul(is0_eps, is0_eps)));
			const vector4double inv_asum = myreciprocal<preclevel>(vec_add(a0, a1));
			const vector4double term0 = vec_msub(F_3, c, b);
            vec_sta(vec_max(min_val_minus, vec_min(max_val_minus, vec_mul(vec_mul(F_1_2, inv_asum), vec_madd(a0, term0, vec_mul(a1, term1))))), 0L, minus);
		}

		//plus. nothing to share with minus here
		{
			const vector4double dme = vec_sub(d, e);
			const vector4double is0p_eps = vec_madd(dme, dme, mywenoeps);
			const vector4double a0p = vec_mul(F_1_3, myreciprocal<preclevel>(vec_mul(is0p_eps, is0p_eps)));
			const vector4double inv_asump = myreciprocal<preclevel>(vec_add(a0p, a1));
			const vector4double term0p = vec_msub(F_3, d, e);
            vec_sta(vec_max(min_val_plus, vec_min(max_val_plus, vec_mul(vec_mul(F_1_2, inv_asump), vec_madd(a0p, term0p, vec_mul(a1, term1))))), 0L, plus);
		}
	}
};

template<>
template<>
inline void WenoQPXUnrolled<WenoQPX_3rdOrder_Minus>::_weno_unrolled<0>(Real * const a, Real * const b, Real * const c,
																	 Real * const d, Real * const e, Real * const out) const{}

template<>
template<>
inline void WenoQPXUnrolled<WenoQPX_3rdOrder_Plus>::_weno_unrolled<0>(Real * const a, Real * const b, Real * const c,
																	Real * const d, Real * const e, Real * const out) const{}
