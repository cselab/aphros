#pragma once

#ifdef _QPX_
#define vec_min(a, b)   vec_sel(a, b, vec_cmpgt(a, b))
#define vec_max(a, b)   vec_sel(a, b, vec_cmplt(a, b))
#endif

struct Weno_QPX_fused
{
	vector4double F_4_10, F_11_10, M_19_10, F_25_10, M_31_10;		/* A */
	vector4double F_4_13, F_5_13;									/* B */
	vector4double F_10_3, F_13_3;									/* C */
	vector4double M_1_2, F_5_2, M_7_2, F_11_2;						/* D */
	vector4double F_2_6;											/* E */
	vector4double WENOEPS4;

	Weno_QPX_fused()
	{
		/* A */
		F_4_10 = vec_splats(4.f/10.f);
		F_11_10 = vec_splats(11.f/10.f);
		M_19_10 = vec_splats(-19.f/10.f);
		F_25_10 = vec_splats(25.f/10.f);
		M_31_10 = vec_splats(-31.f/10.f);

		/* B */
		F_4_13 = vec_splats(4.f/13.f);
		F_5_13 = vec_splats(5.f/13.f);

		/* C */
		F_10_3 = vec_splats(10.f/3.f);
		F_13_3 = vec_splats(13.f/3.f);

		/* D */
		M_1_2 = vec_splats(-1.f/2.f);
		F_5_2 = vec_splats(5.f/2.f);
		M_7_2 = vec_splats(-7.f/2.f);
		F_11_2 = vec_splats(11.f/2.f);

		/* E */
		F_2_6 = vec_splats(2.f/6.f);

		WENOEPS4 = vec_splats(1e-6f);
	}

	void weno_minus_plus_fused_opt1(const vector4double a, const vector4double b, const vector4double c, const vector4double d, const vector4double e, const vector4double f, Real * minus, Real * plus) //133 FLOP
	{
        vector4double min_val_minus = vec_min(a, vec_min(b, vec_min(c, vec_min(d,e))));
        vector4double max_val_minus = vec_max(a, vec_max(b, vec_max(c, vec_max(d,e))));

        vector4double min_val_plus = vec_min(f, vec_min(b, vec_min(c, vec_min(d,e))));
        vector4double max_val_plus = vec_max(f, vec_max(b, vec_max(c, vec_max(d,e))));
		vector4double m_is0, m_is1, m_is2;
		vector4double p_is0, p_is1, p_is2;

		const vector4double cd = vec_mul(c, d);	// 1
		const vector4double c2 = vec_mul(c, c);	// 1
		const vector4double d2 = vec_mul(d, d);	// 1

		{
			const vector4double b2 = vec_mul(b, b);	// 1
			const vector4double bc = vec_mul(b, c);	// 1
			const vector4double bd = vec_mul(b, d);	// 1

			{
				const vector4double t1 = vec_madd(a, F_4_10, vec_madd(b, M_19_10, vec_mul(c, F_11_10)));
				const vector4double t2 = vec_madd(b2, F_25_10, vec_madd(bc, M_31_10, c2));
				m_is0 = vec_madd(a, t1, t2);
			}

			{
				const vector4double t1 = vec_msub(b2, F_4_13, bc);
				const vector4double t2 = vec_madd(bd, F_5_13, c2);
				const vector4double t3 = vec_msub(d2, F_4_13, cd);
				m_is1 = vec_add(t1, vec_add(t2, t3));
			}

			{
				p_is2 = vec_madd(b2, F_4_10, vec_madd(bc, M_19_10, vec_madd(c2, F_25_10, vec_madd(bd, F_11_10, vec_madd(cd, M_31_10, d2)))));
			}
		}

		{
			const vector4double ce = vec_mul(c, e);	// 1
			const vector4double de = vec_mul(d, e);	// 1
			const vector4double e2 = vec_mul(e, e);	// 1

			{
				m_is2 = vec_madd(e2, F_4_10, vec_madd(de, M_19_10, vec_madd(d2, F_25_10, vec_madd(ce, F_11_10, vec_madd(cd, M_31_10, c2)))));
			}

			{
				const vector4double t1 = vec_msub(c2, F_4_13, cd);
				const vector4double t2 = vec_madd(ce, F_5_13, d2);
				const vector4double t3 = vec_msub(e2, F_4_13, de);
				p_is1 = vec_add(t1, vec_add(t2, t3));
			}

			{
				const vector4double t1 = vec_madd(f, F_4_10, vec_madd(e, M_19_10, vec_mul(d, F_11_10)));
				const vector4double t2 = vec_madd(e2, F_25_10, vec_madd(de, M_31_10, d2));
				p_is0 = vec_madd(f, t1, t2);
			}
		}
		// 67

		const vector4double m_is0plus = vec_madd(F_10_3, m_is0, WENOEPS4);
		const vector4double m_is1plus = vec_madd(F_13_3, m_is1, WENOEPS4);
		const vector4double m_is2plus = vec_madd(F_10_3, m_is2, WENOEPS4);

		const vector4double p_is0plus = vec_madd(F_10_3, p_is0, WENOEPS4);
		const vector4double p_is1plus = vec_madd(F_13_3, p_is1, WENOEPS4);
		const vector4double p_is2plus = vec_madd(F_10_3, p_is2, WENOEPS4);
		// 12

		const vector4double m_alpha0 = mydivision<preclevel>(vec_splats(0.1f), vec_mul(m_is0plus, m_is0plus));
		const vector4double m_alpha1 = mydivision<preclevel>(vec_splats(0.6f), vec_mul(m_is1plus, m_is1plus));
		const vector4double m_alpha2 = mydivision<preclevel>(vec_splats(0.3f), vec_mul(m_is2plus, m_is2plus));
		const vector4double m_inv_alpha = myreciprocal<preclevel>(vec_add(m_alpha0, vec_add(m_alpha1, m_alpha2)));

		const vector4double m_omega0 = vec_mul(m_alpha0, m_inv_alpha);
		const vector4double m_omega1 = vec_mul(m_alpha1, m_inv_alpha);
		const vector4double m_omega2 = vec_sub(vec_splats(1), vec_add(m_omega0, m_omega1));

		const vector4double p_alpha0 = mydivision<preclevel>(vec_splats(0.1f), vec_mul(p_is0plus, p_is0plus));
		const vector4double p_alpha1 = mydivision<preclevel>(vec_splats(0.6f), vec_mul(p_is1plus, p_is1plus));
		const vector4double p_alpha2 = mydivision<preclevel>(vec_splats(0.3f), vec_mul(p_is2plus, p_is2plus));
		const vector4double p_inv_alpha = myreciprocal<preclevel>(vec_add(p_alpha0, vec_add(p_alpha1, p_alpha2)));

		const vector4double p_omega0 = vec_mul(p_alpha0, p_inv_alpha);
		const vector4double p_omega1 = vec_mul(p_alpha1, p_inv_alpha);
		const vector4double p_omega2 = vec_sub(vec_splats(1), vec_add(p_omega0, p_omega1));
		// 26

		{
			const vector4double m1_p2 = vec_madd(b, M_1_2, vec_madd(c, F_5_2, d)); //M_1_2*b + F_5_2*c + d;	// 4
			const vector4double m2_p1 = vec_madd(e, M_1_2, vec_madd(d, F_5_2, c));//c + F_5_2*d + M_1_2*e;		// 4
			const vector4double m0 = vec_madd(b, M_7_2, vec_madd(c, F_11_2, a)); //a +M_7_2*b + F_11_2*c
			const vector4double p0 = vec_madd(e, M_7_2, vec_madd(d, F_11_2, f)); //f +M_7_2*e + F_11_2*d

			const vector4double minus_tmp = vec_madd(m_omega0, m0, vec_madd(m_omega1, m1_p2, vec_mul(m_omega2, m2_p1)));
			const vector4double plus_tmp  = vec_madd(p_omega0, p0, vec_madd(p_omega1, m2_p1, vec_mul(p_omega2, m1_p2)));

                        vec_sta(vec_max(min_val_minus, vec_min(max_val_minus, vec_mul(minus_tmp, F_2_6))), 0L, minus);
                        vec_sta(vec_max(min_val_plus, vec_min(max_val_plus, vec_mul(plus_tmp, F_2_6))), 0L, plus);
		}
		// 28
	}

	void weno_minus_plus_fused_opt2(const vector4double a, const vector4double b, const vector4double c, const vector4double d, const vector4double e, const vector4double f, Real *minus, Real *plus) //127 FLOP
	{
        vector4double min_val_minus = vec_min(a, vec_min(b, vec_min(c, vec_min(d,e))));
        vector4double max_val_minus = vec_max(a, vec_max(b, vec_max(c, vec_max(d,e))));

        vector4double min_val_plus = vec_min(f, vec_min(b, vec_min(c, vec_min(d,e))));
        vector4double max_val_plus = vec_max(f, vec_max(b, vec_max(c, vec_max(d,e))));

		///
		vector4double m_is0, m_is1, m_is2;
		vector4double p_is0, p_is1, p_is2;

		const vector4double cd = vec_mul(c, d);	// 1
		const vector4double c2 = vec_mul(c, c);	// 1
		const vector4double d2 = vec_mul(d, d);	// 1

		{
			const vector4double b2 = vec_mul(b, b);	// 1
			const vector4double bc = vec_mul(b, c);	// 1
			const vector4double bd = vec_mul(b, d);	// 1

			{
				const vector4double t1 = vec_madd(a, F_4_10, vec_madd(b, M_19_10, vec_mul(c, F_11_10)));
				const vector4double t2 = vec_madd(b2, F_25_10, vec_madd(bc, M_31_10, c2));
				m_is0 = vec_madd(a, t1, t2);
			}

			{
				const vector4double t1 = vec_msub(b2, F_4_13, bc);
				const vector4double t2 = vec_madd(bd, F_5_13, c2);
				const vector4double t3 = vec_msub(d2, F_4_13, cd);
				m_is1 = vec_add(t1, vec_add(t2, t3));
			}

			{
				p_is2 = vec_madd(b2, F_4_10, vec_madd(bc, M_19_10, vec_madd(c2, F_25_10, vec_madd(bd, F_11_10, vec_madd(cd, M_31_10, d2)))));
			}
		}

		{
			const vector4double ce = vec_mul(c, e);	// 1
			const vector4double de = vec_mul(d, e);	// 1
			const vector4double e2 = vec_mul(e, e);	// 1

			{
				m_is2 = vec_madd(e2, F_4_10, vec_madd(de, M_19_10, vec_madd(d2, F_25_10, vec_madd(ce, F_11_10, vec_madd(cd, M_31_10, c2)))));
			}

			{
				const vector4double t1 = vec_msub(c2, F_4_13, cd);
				const vector4double t2 = vec_madd(ce, F_5_13, d2);
				const vector4double t3 = vec_msub(e2, F_4_13, de);
				p_is1 = vec_add(t1, vec_add(t2, t3));
			}

			{
				const vector4double t1 = vec_madd(f, F_4_10, vec_madd(e, M_19_10, vec_mul(d, F_11_10)));
				const vector4double t2 = vec_madd(e2, F_25_10, vec_madd(de, M_31_10, d2));
				p_is0 = vec_madd(f, t1, t2);
			}
		}
		// 67

		const vector4double m_is0plus = vec_madd(F_10_3, m_is0, WENOEPS4);
		const vector4double m_is1plus = vec_madd(F_13_3, m_is1, WENOEPS4);
		const vector4double m_is2plus = vec_madd(F_10_3, m_is2, WENOEPS4);

		const vector4double p_is0plus = vec_madd(F_10_3, p_is0, WENOEPS4);
		const vector4double p_is1plus = vec_madd(F_13_3, p_is1, WENOEPS4);
		const vector4double p_is2plus = vec_madd(F_10_3, p_is2, WENOEPS4);

		const vector4double m_alpha0 = mydivision<preclevel>(vec_splats(0.1f), vec_mul(m_is0plus, m_is0plus));
		const vector4double m_alpha1 = mydivision<preclevel>(vec_splats(0.6f), vec_mul(m_is1plus, m_is1plus));
		const vector4double m_alpha2 = mydivision<preclevel>(vec_splats(0.3f), vec_mul(m_is2plus, m_is2plus));
		const vector4double m_inv_alpha = myreciprocal<preclevel>(vec_add(m_alpha0, vec_add(m_alpha1, m_alpha2)));

		const vector4double p_alpha0 = mydivision<preclevel>(vec_splats(0.1f), vec_mul(p_is0plus, p_is0plus));
		const vector4double p_alpha1 = mydivision<preclevel>(vec_splats(0.6f), vec_mul(p_is1plus, p_is1plus));
		const vector4double p_alpha2 = mydivision<preclevel>(vec_splats(0.3f), vec_mul(p_is2plus, p_is2plus));
		const vector4double p_inv_alpha = myreciprocal<preclevel>(vec_add(p_alpha0, vec_add(p_alpha1, p_alpha2)));

		{
			const vector4double m1_p2 = vec_madd(b, M_1_2, vec_madd(c, F_5_2, d));
			const vector4double m2_p1 = vec_madd(e, M_1_2, vec_madd(d, F_5_2, c));
			const vector4double m0 = vec_madd(b, M_7_2, vec_madd(c, F_11_2, a));
			const vector4double p0 = vec_madd(e, M_7_2, vec_madd(d, F_11_2, f));

			const vector4double minus_tmp = vec_madd(m_alpha0, m0, vec_madd(m_alpha1, m1_p2, vec_mul(m_alpha2, m2_p1)));
			const vector4double plus_tmp  = vec_madd(p_alpha0, p0, vec_madd(p_alpha1, m2_p1, vec_mul(p_alpha2, m1_p2)));

                        vec_sta(vec_max(min_val_minus, vec_min(max_val_minus, vec_mul(minus_tmp, vec_mul(F_2_6, m_inv_alpha)))), 0L, minus);
                        vec_sta(vec_max(min_val_plus, vec_min(max_val_plus, vec_mul(plus_tmp , vec_mul(F_2_6, p_inv_alpha)))), 0L, plus);
		}
	}
};
