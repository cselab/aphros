/*
 *  *  WenoQPX.h
 *   *  TryThis
 *    *
 *     *  Created by Diego Rossinelli on 1/24/13.
 *      *  Copyright 2013 ETH Zurich. All rights reserved.
 *       *
 *        */

#pragma once
#ifdef _QPXEMU_
#include "QPXEMU.h"
#endif

struct WenoQPX_BaseFunctor
{
protected:
	vector4double one, F_4_3, F_5_3, F_10_3, F_11_3,  F_13_3,  F_19_3, F_25_3, F_31_3; 
	vector4double mywenoeps;
	vector4double M_1_6, F_1_3, F_5_6, M_7_6, F_11_6; 
	//vector4double F_1_10, F_6_10, F_3_10;
	
	inline vector4double _compute_is0(const vector4double a, const vector4double b, const vector4double c) const
	{
		const vector4double x = vec_madd(a, F_4_3, vec_msub(c, F_11_3, vec_mul(b, F_19_3)));
		const vector4double y = vec_msub(b, F_25_3, vec_mul(c, F_31_3));
		const vector4double z = vec_mul(vec_mul(c, c), F_10_3);
		return vec_madd(a, x, vec_madd(b, y, z));
	}
	
	inline vector4double _compute_is1(const vector4double b, const vector4double c, const vector4double d) const
	{
		const vector4double x = vec_madd(b, F_4_3, vec_msub(d, F_5_3, vec_mul(c, F_13_3)));
		const vector4double y = vec_msub(c, F_13_3, vec_mul(d, F_13_3));
		const vector4double z = vec_mul(vec_mul(d, d), F_4_3);
		return vec_madd(b, x, vec_madd(c, y, z));
	}
	
	inline vector4double _compute_is2(const vector4double c, const vector4double d, const vector4double e) const
	{
		const vector4double x = vec_madd(c, F_10_3, vec_msub(e, F_11_3, vec_mul(d, F_31_3)));
		const vector4double y = vec_msub(d, F_25_3, vec_mul(e, F_19_3));
		const vector4double z = vec_mul(vec_mul(e, e), F_4_3);
		return vec_madd(c, x, vec_madd(d, y, z));
	}
	
	WenoQPX_BaseFunctor()
	{
		one = vec_splats(1.); 
		F_4_3 = vec_splats(4./3.); 
		F_5_3 = vec_splats(5./3.); 
		F_10_3 = vec_splats(10./3.); 
		F_11_3 = vec_splats(11./3.); 
		F_13_3 = vec_splats(13./3.); 
		F_19_3 = vec_splats(19./3.); 
		F_25_3 = vec_splats(25./3.); 
		F_31_3 = vec_splats(31./3.); 
		mywenoeps = vec_splats(WENOEPS); 
		
		M_1_6 = vec_splats(-1./6.); 
		F_1_3 = vec_splats(1./3.); 
		F_5_6 = vec_splats(5./6.); 
		M_7_6 = vec_splats(-7./6.); 
		F_11_6 = vec_splats(11./6.); 
		
	//	F_1_10 = vec_splats(1./10.); 
	//	F_6_10 = vec_splats(6./10.); 
	//	F_3_10 = vec_splats(3./10.);			
	}	
};

class WenoQPX_MinusFunctor : protected WenoQPX_BaseFunctor
{
public:
	inline vector4double _weno(const vector4double a, const vector4double b, const vector4double c, 
							   const vector4double d, const vector4double e) const 
	{		
		const vector4double is0 = _compute_is0(a, b, c);		
		const vector4double is1 = _compute_is1(b, c, d);
		const vector4double is2 = _compute_is2(c, d, e);
		
		const vector4double is0plus = vec_add(is0, mywenoeps);
		const vector4double is1plus = vec_add(is1, mywenoeps);
		const vector4double is2plus = vec_add(is2, mywenoeps);
		
		const vector4double alpha0 = myreciprocal<preclevel>(vec_mul(is0plus, is0plus));		
		const vector4double alpha1 = mydivision<preclevel>(vec_splats(6), vec_mul(is1plus, is1plus));
		const vector4double alpha2 = mydivision<preclevel>(vec_splats(3), vec_mul(is2plus, is2plus));
		const vector4double inv_alpha = myreciprocal<preclevel>(vec_add(alpha0, vec_add(alpha1, alpha2))); 
		
		const vector4double omega0 = vec_mul(alpha0, inv_alpha);
		const vector4double omega1 = vec_mul(alpha1, inv_alpha);
		const vector4double omega2 = vec_sub(one, vec_add(omega0, omega1));		
		
		const vector4double t0 = vec_madd(a, F_1_3, vec_madd(b, M_7_6,  vec_mul(c, F_11_6)));
		const vector4double t1 = vec_madd(b, M_1_6, vec_madd(c, F_5_6,  vec_mul(d, F_1_3)));
		const vector4double t2 = vec_madd(c, F_1_3, vec_madd(d, F_5_6,  vec_mul(e, M_1_6)));
		
		return vec_madd(omega0, t0, vec_madd(omega1, t1, vec_mul(omega2, t2)));	
	}
};

class WenoQPX_PlusFunctor : protected WenoQPX_BaseFunctor
{
public:
	inline vector4double _weno(const vector4double b, const vector4double c, const vector4double d, 
							   const vector4double e, const vector4double f) const 
	{
		const vector4double is0 = _compute_is2(d, e, f);		
		const vector4double is1 = _compute_is1(c, d, e);
		const vector4double is2 = _compute_is0(b, c, d);
		
		const vector4double is0plus = vec_add(is0, mywenoeps);
		const vector4double is1plus = vec_add(is1, mywenoeps);
		const vector4double is2plus = vec_add(is2, mywenoeps);
		
		const vector4double alpha0 = myreciprocal<preclevel>(vec_mul(is0plus, is0plus));		
		const vector4double alpha1 = mydivision<preclevel>(vec_splats(6), vec_mul(is1plus, is1plus));
		const vector4double alpha2 = mydivision<preclevel>(vec_splats(3), vec_mul(is2plus, is2plus));
		const vector4double inv_alpha = myreciprocal<preclevel>(vec_add(alpha0, vec_add(alpha1, alpha2))); 
		
		const vector4double omega0 = vec_mul(alpha0, inv_alpha);
		const vector4double omega1 = vec_mul(alpha1, inv_alpha);
		const vector4double omega2 = vec_sub(one, vec_add(omega0, omega1));		
		
		const vector4double t0 = vec_madd(f, F_1_3, vec_madd(e, M_7_6,  vec_mul(d, F_11_6)));
		const vector4double t1 = vec_madd(e, M_1_6, vec_madd(d, F_5_6,  vec_mul(c, F_1_3)));
		const vector4double t2 = vec_madd(d, F_1_3, vec_madd(c, F_5_6,  vec_mul(b, M_1_6)));
		
		return vec_madd(omega0, t0, vec_madd(omega1, t1, vec_mul(omega2, t2)));	
	}
};

struct WenoQPX_PanosFunctor
{
protected:
	
	vector4double one, F_4_3, F_5_3, F_10_3, F_11_3,  F_13_3,  F_19_3, F_25_3, F_31_3; 
	vector4double mywenoeps;
	//screw this vector4double M_1_6, F_1_3, F_5_6, M_7_6, F_11_6; 
	//screw this vector4double F_6_10, F_3_10;
	vector4double M_1_2, F_1_3, F_5_2, M_7_2, F_11_2;
	
	inline vector4double _compute_is0(const vector4double zpanos, const vector4double a, const vector4double b, const vector4double c) const
	{
		const vector4double x = vec_madd(a, F_4_3, vec_msub(c, F_11_3, vec_mul(b, F_19_3)));
		const vector4double y = vec_msub(b, F_25_3, vec_mul(c, F_31_3));
		//screw this const vector4double z = vec_mul(vec_mul(c, c), F_10_3);
		return vec_madd(a, x, vec_madd(b, y, zpanos));
	}
	
	inline vector4double _compute_is1(const vector4double b, const vector4double c, const vector4double d) const
	{
		const vector4double x = vec_madd(b, F_4_3, vec_msub(d, F_5_3, vec_mul(c, F_13_3)));
		const vector4double y = vec_msub(c, F_13_3, vec_mul(d, F_13_3));
		const vector4double z = vec_mul(vec_mul(d, d), F_4_3);
		return vec_madd(b, x, vec_madd(c, y, z));
	}
	
	inline vector4double _compute_is2(const vector4double zpanos, const vector4double c, const vector4double d, const vector4double e) const
	{
		//screw that const vector4double x = vec_madd(c, F_10_3, vec_msub(e, F_11_3, vec_mul(d, F_31_3)));
		const vector4double x = vec_msub(e, F_11_3, vec_mul(d, F_31_3));
		const vector4double y = vec_msub(d, F_25_3, vec_mul(e, F_19_3));
		const vector4double z = vec_mul(vec_mul(e, e), F_4_3);
		//screw this return vec_madd(c, x, vec_madd(d, y, z));
		return vec_add(zpanos, vec_madd(c, x, vec_madd(d, y, z)));
	}
	
	WenoQPX_PanosFunctor()
	{
		one = vec_splats(1.); 
		F_4_3 = vec_splats(4./3.); 
		F_5_3 = vec_splats(5./3.); 
		F_10_3 = vec_splats(10./3.); 
		F_11_3 = vec_splats(11./3.); 
		F_13_3 = vec_splats(13./3.); 
		F_19_3 = vec_splats(19./3.); 
		F_25_3 = vec_splats(25./3.); 
		F_31_3 = vec_splats(31./3.); 
		mywenoeps = vec_splats(WENOEPS); 
		
		M_1_2 = vec_splats(-0.5f);
		F_1_3 = vec_splats(1./3.);
		F_5_2 = vec_splats(2.5f);
		M_7_2 = vec_splats(-3.5f); 
		F_11_2 = vec_splats(5.5f);
	//	M_1_6 = vec_splats(-1./6.); 
//		F_1_3 = vec_splats(1./3.); 
//		F_5_6 = vec_splats(5./6.); 
//		M_7_6 = vec_splats(-7./6.); 
//		F_11_6 = vec_splats(11./6.); 
//		
	//	F_1_10 = vec_splats(1./10.); 
	//	F_6_10 = vec_splats(6./10.); 
	//	F_3_10 = vec_splats(3./10.);			
	}	
};

class WenoQPX_MinusPanos : protected WenoQPX_PanosFunctor
{
public:
	
	inline vector4double _weno(const vector4double a, const vector4double b, const vector4double c, 
							   const vector4double d, const vector4double e) const 
	{		
		const vector4double zpanos = vec_mul(c, vec_mul(c, F_10_3));
		const vector4double is0 = _compute_is0(zpanos, a, b, c);		
		const vector4double is1 = _compute_is1(b, c, d);
		const vector4double is2 = _compute_is2(zpanos, c, d, e);
		
		const vector4double is0plus = vec_add(is0, mywenoeps);
		const vector4double is1plus = vec_add(is1, mywenoeps);
		const vector4double is2plus = vec_add(is2, mywenoeps);
		
		const vector4double alpha0 = myreciprocal<preclevel>(vec_mul(is0plus, is0plus));		
		const vector4double alpha1 = mydivision<preclevel>(vec_splats(6), vec_mul(is1plus, is1plus));
		const vector4double alpha2 = mydivision<preclevel>(vec_splats(3), vec_mul(is2plus, is2plus));
		const vector4double inv_alpha = myreciprocal<preclevel>(vec_add(alpha0, vec_add(alpha1, alpha2))); 
		
		//const vector4double omega0 = vec_mul(alpha0, inv_alpha);
		//const vector4double omega1 = vec_mul(alpha1, inv_alpha);
		//const vector4double omega2 = vec_sub(one, vec_add(omega0, omega1));		
		
		const vector4double t0 = vec_madd(b, vec_splats(-3.5f),  vec_madd(c, vec_splats(5.5f), a));
		const vector4double t1 = vec_madd(b, vec_splats(-0.5f),  vec_madd(c, vec_splats(2.5f), d));
		const vector4double t2 = vec_madd(d, vec_splats(2.5f),  vec_madd(e, vec_splats(-0.5f), c));
		
		const vector4double panos2 = vec_madd(alpha0, t0, vec_madd(alpha1, t1, vec_mul(alpha2, t2)));
		return vec_mul(vec_mul(vec_splats(1.f/3.f), inv_alpha), panos2);
	}
};

class WenoQPX_PlusPanos : protected WenoQPX_PanosFunctor
{
public:
	
	inline vector4double _weno(const vector4double b, const vector4double c, const vector4double d, 
							   const vector4double e, const vector4double f) const 
	{
		const vector4double zpanos = vec_mul(d, vec_mul(d, F_10_3));
		const vector4double is0 = _compute_is2(zpanos, d, e, f);		
		const vector4double is1 = _compute_is1(c, d, e);
		const vector4double is2 = _compute_is0(zpanos, b, c, d);
		
		const vector4double is0plus = vec_add(is0, mywenoeps);
		const vector4double is1plus = vec_add(is1, mywenoeps);
		const vector4double is2plus = vec_add(is2, mywenoeps);
		
		const vector4double alpha0 = myreciprocal<preclevel>(vec_mul(is0plus, is0plus));		
		const vector4double alpha1 = mydivision<preclevel>(vec_splats(6), vec_mul(is1plus, is1plus));
		const vector4double alpha2 = mydivision<preclevel>(vec_splats(3), vec_mul(is2plus, is2plus));
		const vector4double inv_alpha = myreciprocal<preclevel>(vec_add(alpha0, vec_add(alpha1, alpha2))); 
		
		//const vector4double omega0 = vec_mul(alpha0, inv_alpha);
		//const vector4double omega1 = vec_mul(alpha1, inv_alpha);
		//const vector4double omega2 = vec_sub(one, vec_add(omega0, omega1));		
		
		const vector4double t0 = vec_madd(d, vec_splats(5.5f), vec_madd(e, vec_splats(-3.5f), f));
		const vector4double t1 = vec_madd(e, vec_splats(-0.5f), vec_madd(d, vec_splats(2.5f), c));
		const vector4double t2 = vec_madd(c, vec_splats(2.5f), vec_madd(b, vec_splats(-0.5f), d));
		
		//return vec_madd(omega0, t0, vec_madd(omega1, t1, vec_mul(omega2, t2)));	
		const vector4double panos2 = vec_madd(alpha0, t0, vec_madd(alpha1, t1, vec_mul(alpha2, t2)));
		return vec_mul(vec_mul(vec_splats(1.f/3.f), inv_alpha), panos2);
	}
};

template<typename MyWenoFunctor>
class WenoQPX : protected MyWenoFunctor
{
public:
	
	void compute(Real * const a, Real * const b, Real * const c, 
				 Real * const d, Real * const e, Real * const out,
				 const int NENTRIES) const 
	{
		enum { size_of_real = sizeof(Real) } ;
		
#pragma omp for schedule(runtime) 
		for(int i=0; i<NENTRIES; i += 4)
		{
			const unsigned long offset = size_of_real * i;
			
			const vector4double mya = vec_lda(offset, a);
			const vector4double myb = vec_lda(offset, b);
			const vector4double myc = vec_lda(offset, c);
			const vector4double myd = vec_lda(offset, d);
			const vector4double mye = vec_lda(offset, e);
			
			const vector4double result = this->_weno(mya, myb, myc, myd, mye);
			
			vec_sta(result, offset, out);
		}
	}
	
	template<int NENTRIES>
	void stride(Real * const a, Real * const b, Real * const c, 
				Real * const d, Real * const e, Real * const out) const 
	{
		enum { size_of_real = sizeof(Real) } ;
				
#pragma unroll(4)
		for(int i=0; i<NENTRIES; i += 4)
		{
			const unsigned long offset = size_of_real * i;
			
			const vector4double mya = vec_lda(offset, a);
			const vector4double myb = vec_lda(offset, b);
			const vector4double myc = vec_lda(offset, c);
			const vector4double myd = vec_lda(offset, d);
			const vector4double mye = vec_lda(offset, e);
			
			const vector4double result = this->_weno(mya, myb, myc, myd, mye);
			
			vec_sta(result, offset, out);
		}
		//#endif
	}
	
	void operator()(const Real * const a, const Real * const b, const Real * const c, 
					const Real * const d, const Real * const e, Real * const out,
					const int NENTRIES) const 
	{
		compute(const_cast<Real * const>(a),
				const_cast<Real * const>(b),
				const_cast<Real * const>(c),
				const_cast<Real * const>(d),
				const_cast<Real * const>(e),
				out, NENTRIES);
	}
};

template<typename MyWenoFunctor>
class WenoQPXUnrolled : protected WenoQPX<MyWenoFunctor>
{
#ifndef _WENO_UNROLLS_
#define _WENO_UNROLLS_ 4
#endif
	
	enum { UNROLLS = _WENO_UNROLLS_ };
	
public:
	template<int count>
	inline void _weno_unrolled(Real * const a, Real * const b, Real * const c, 
							   Real * const d, Real * const e, Real * const out) const
	{
		enum 
		{ 
			size_of_real = sizeof(Real),
			offset = size_of_real * 4 * (UNROLLS - count)
		} ;
		
		const vector4double result = 
		WenoQPX<MyWenoFunctor>::_weno(vec_lda(offset, a), vec_lda(offset, b), 
									  vec_lda(offset, c), vec_lda(offset, d), 
									  vec_lda(offset, e));
		
		vec_sta(result, offset, out);
		
		WenoQPXUnrolled::template _weno_unrolled<count - 1>(a, b, c, d, e, out);
	}
	
	WenoQPXUnrolled() : WenoQPX<MyWenoFunctor>() {}
	
	template<int NENTRIES >
	void stride(Real * const a, Real * const b, Real * const c, 
				Real * const d, Real * const e, Real * const out) const 
	{
		enum { 
			UNROLLED_ITERS = 4 * UNROLLS ,
			nice = UNROLLED_ITERS * (NENTRIES / UNROLLED_ITERS)
		};
		
		/* #pragma omp for schedule(runtime) */
		for(int i = 0; i < nice; i += UNROLLED_ITERS)
			_weno_unrolled<UNROLLS>(a + i, b + i, c + i, d + i, e + i, out + i);
		
		WenoQPX<MyWenoFunctor>::compute(a + nice, b + nice, c + nice, d + nice, e + nice, 
										out + nice, NENTRIES - nice);	
	}
	
	
	void compute(Real * const a, Real * const b, Real * const c, 
				 Real * const d, Real * const e, Real * const out,
				 const int NENTRIES) const 
	{
		enum { UNROLLED_ITERS = 4 * UNROLLS };
		const int nice = UNROLLED_ITERS * (NENTRIES / UNROLLED_ITERS);
		
#pragma omp for schedule(runtime)
		for(int i = 0; i < nice; i += UNROLLED_ITERS)
			_weno_unrolled<UNROLLS>(a + i, b + i, c + i, d + i, e + i, out + i);
		
		WenoQPX<MyWenoFunctor>::compute(a + nice, b + nice, c + nice, d + nice, e + nice, 
										out + nice, NENTRIES - nice);	
	}
	
	void operator()(const Real * const a, const Real * const b, const Real * const c, 
					const Real * const d, const Real * const e, Real * const out,
					const int NENTRIES) const 
	{
		compute(const_cast<Real * const>(a),
				const_cast<Real * const>(b),
				const_cast<Real * const>(c),
				const_cast<Real * const>(d),
				const_cast<Real * const>(e),
				out, NENTRIES);
	}
};

template<>
template<>
inline void WenoQPXUnrolled<WenoQPX_MinusFunctor>::_weno_unrolled<0>(Real * const a, Real * const b, Real * const c, 
																	 Real * const d, Real * const e, Real * const out) const{}

template<>
template<>
inline void WenoQPXUnrolled<WenoQPX_PlusFunctor>::_weno_unrolled<0>(Real * const a, Real * const b, Real * const c, 
																	Real * const d, Real * const e, Real * const out) const{}

