/*
 *  *  Weno_CPP.h
 *   *  TryThis
 *    *
 *     *  Created by Diego Rossinelli on 1/24/13.
 *      *  Copyright 2013 ETH Zurich. All rights reserved.
 *       *
 *        */

#pragma once

#ifdef __xlC__
extern
#ifdef __cplusplus
"builtin"
#endif
void __alignx (int n, const void *addr);
#endif

struct Weno_FunctorMinus
{
	inline Real _weno(const Real a, const Real b, const Real c, const Real d, const Real e) const //82 FLOP
	{	
#ifdef _WENO3_
		const Real is0 = (c-b)*(c-b);
		const Real is1 = (d-c)*(d-c);
		
		const Real alpha0 = 1./(3.*(is0+WENOEPS)*(is0+WENOEPS));
		const Real alpha1 = 2./(3.*(is1+WENOEPS)*(is1+WENOEPS));
		
		const Real omega0=alpha0/(alpha0+alpha1);
		const Real omega1=1.-omega0;
		
		return omega0*(1.5*c-.5*b) + omega1*(.5*c+.5*d);
#else
		const Real is0 = a*(a*(Real)(4./3.)  - b*(Real)(19./3.)  + c*(Real)(11./3.)) + b*(b*(Real)(25./3.)  - c*(Real)(31./3.)) + c*c*(Real)(10./3.);		
		const Real is1 = b*(b*(Real)(4./3.)  - c*(Real)(13./3.)  + d*(Real)(5./3.))  + c*(c*(Real)(13./3.)  - d*(Real)(13./3.)) + d*d*(Real)(4./3.);
		const Real is2 = c*(c*(Real)(10./3.) - d*(Real)(31./3.)  + e*(Real)(11./3.)) + d*(d*(Real)(25./3.)  - e*(Real)(19./3.)) + e*e*(Real)(4./3.);
		
		const Real is0plus = is0 + (Real)WENOEPS;
		const Real is1plus = is1 + (Real)WENOEPS;
		const Real is2plus = is2 + (Real)WENOEPS;
		
		const Real alpha0 = (Real)(0.1)*((Real)1/(is0plus*is0plus));
		const Real alpha1 = (Real)(0.6)*((Real)1/(is1plus*is1plus));
		const Real alpha2 = (Real)(0.3)*((Real)1/(is2plus*is2plus));
		const Real alphasum = alpha0+alpha1+alpha2;
		const Real inv_alpha = ((Real)1)/alphasum;
		
		const Real omega0 = alpha0 * inv_alpha;
		const Real omega1 = alpha1 * inv_alpha;
		const Real omega2 = 1-omega0-omega1;
		
		return omega0*((Real)(1.0/3.)*a-(Real)(7./6.)*b+(Real)(11./6.)*c) + omega1*(-(Real)(1./6.)*b+(Real)(5./6.)*c+(Real)(1./3.)*d) + omega2*((Real)(1./3.)*c+(Real)(5./6.)*d-(Real)(1./6.)*e);
#endif
	}
};

struct Weno_FunctorPlus
{
	inline Real _weno(const Real b, const Real c, const Real d, const Real e, const Real f) const//82 FLOP
	{	
#ifdef _WENO3_
		const Real is0 = (d-e)*(d-e);
		const Real is1 = (d-c)*(d-c);
		
		const Real alpha0 = (1./3.)/((is0+WENOEPS)*(is0+WENOEPS));
		const Real alpha1 = (2./3.)/((is1+WENOEPS)*(is1+WENOEPS));
		
		const Real omega0 = alpha0/(alpha0+alpha1);
		const Real omega1 = 1.-omega0;
		
		return omega0*(1.5*d-.5*e) + omega1*(.5*d+.5*c);
#else
		const Real is0 = d*(d*(Real)(10./3.)- e*(Real)(31./3.) + f*(Real)(11./3.)) + e*(e*(Real)(25./3.) - f*(Real)(19./3.)) +	f*f*(Real)(4./3.);
		const Real is1 = c*(c*(Real)(4./3.) - d*(Real)(13./3.) + e*(Real)(5./3.)) + d*(d*(Real)(13./3.)  - e*(Real)(13./3.)) +	e*e*(Real)(4./3.);
		const Real is2 = b*(b*(Real)(4./3.) - c*(Real)(19./3.) + d*(Real)(11./3.)) + c*(c*(Real)(25./3.) - d*(Real)(31./3.)) +	d*d*(Real)(10./3.);
		
		const Real is0plus = is0 + (Real)WENOEPS;
		const Real is1plus = is1 + (Real)WENOEPS;
		const Real is2plus = is2 + (Real)WENOEPS;
		
		const Real alpha0 = (Real)(0.1)*(((Real)1)/(is0plus*is0plus));
		const Real alpha1 = (Real)(0.6)*(((Real)1)/(is1plus*is1plus));
		const Real alpha2 = (Real)(0.3)*(((Real)1)/(is2plus*is2plus));
		const Real alphasum = alpha0+alpha1+alpha2;
		
		const Real omega0=alpha0 * (((Real)1)/alphasum);
		const Real omega1=alpha1 * (((Real)1)/alphasum);
		const Real omega2= 1-omega0-omega1;
		
		return omega0*((Real)(1./3.)*f-(Real)(7./6.)*e+(Real)(11./6.)*d) + omega1*(-(Real)(1./6.)*e+(Real)(5./6.)*d+(Real)(1./3.)*c) + omega2*((Real)(1./3.)*d+(Real)(5./6.)*c-(Real)(1./6.)*b);
#endif
	}
	
};

template<typename MyWenoFunctor>
struct WenoReference : MyWenoFunctor
{
	void operator()(const Real * const a, const Real * const b, const Real * const c, 
			  const Real * const d, const Real * const e, Real * const out,
			  const int NENTRIES) const 
	{
#ifdef __xlC__
		__alignx(32, a);
		__alignx(32, b);
		__alignx(32, c);
		__alignx(32, d);
		__alignx(32, e);
		__alignx(32, out);
		
#endif
		
#pragma omp for schedule(runtime)  
		for(int i=0; i<NENTRIES; ++i)
			out[i] = MyWenoFunctor::_weno(a[i], b[i], c[i], d[i], e[i]);
	}	
};

template<typename MyWenoFunctor>
class WenoReferenceUnrolled : protected WenoReference<MyWenoFunctor>
{
#ifndef _WENO_UNROLLS_
#define _WENO_UNROLLS_ 4
#endif
	
	enum 
	{ 
		UNROLLS = _WENO_UNROLLS_,
		UNROLLED_ITERS = 4 * UNROLLS 
	};
	
public:
	template<int count>
	inline void _weno_unrolled(const Real * const a, const Real * const b, const Real * const c, 
							   const Real * const d, const Real * const e, Real * const out) const
	{
		enum { offset = (UNROLLED_ITERS - count) } ;
		
		out[offset] = MyWenoFunctor::_weno(a[offset], b[offset], c[offset], d[offset], e[offset]);
		
		WenoReferenceUnrolled::template _weno_unrolled<count - 1>(a, b, c, d, e, out);
	}
	
	WenoReferenceUnrolled() : WenoReference<MyWenoFunctor>() {}
	
	void operator()(const Real * const a, const Real * const b, const Real * const c, 
				 const Real * const d, const Real * const e, Real * const out,
				 const int NENTRIES) const 
	{
		enum { };
		const int nice = UNROLLED_ITERS * (NENTRIES / UNROLLED_ITERS);
		
#pragma omp for schedule(runtime)
		for(int i = 0; i < nice; i += UNROLLED_ITERS)
			_weno_unrolled<UNROLLED_ITERS>(a + i, b + i, c + i, d + i, e + i, out + i);
		
		WenoReference<MyWenoFunctor>::operator()(a + nice, b + nice, c + nice, d + nice, e + nice, 
										out + nice, NENTRIES - nice);	
	}
};

template<>
template<>
inline void WenoReferenceUnrolled<Weno_FunctorMinus>::_weno_unrolled<0>(const Real * const a, const Real * const b, const Real * const c, 
																		const Real * const d, const Real * const e, Real * const out) const{}

template<>
template<>
inline void WenoReferenceUnrolled<Weno_FunctorPlus>::_weno_unrolled<0>(const Real * const a, const Real * const b, const Real * const c, 
																	   const Real * const d, const Real * const e, Real * const out) const{}


