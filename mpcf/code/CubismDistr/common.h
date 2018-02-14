/*
 *  TwoPhase.h
 *
 *
 *  Created by Diego Rossinelli on 5/14/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

#pragma once
#ifdef _QPXEMU_
#include "QPXEMU.h"
#endif

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif

#ifndef _PREC_LEVEL_
static const int preclevel = 0;
#else
static const int preclevel = _PREC_LEVEL_;
#endif

#include "SOA2D.h"

/*working dataset types */
/*typedef SOA2D<-3, _BLOCKSIZE_+3, -3, _BLOCKSIZE_+3> InputSOA;
typedef RingSOA2D<-3, _BLOCKSIZE_+3, -3, _BLOCKSIZE_+3, 6> RingInputSOA;
typedef SOA2D<0,_BLOCKSIZE_+1, 0,_BLOCKSIZE_> TempSOA;
typedef RingSOA2D<0, _BLOCKSIZE_+1, 0,_BLOCKSIZE_, 2> RingTempSOA;
typedef RingSOA2D<0,_BLOCKSIZE_+1, 0, _BLOCKSIZE_, 3> RingTempSOA3;
*/

//C++ related functions
template<typename X> inline X mysqrt(X x){ abort(); return sqrt(x);}
template<>  inline float mysqrt<float>(float x){ return sqrtf(x);}
template<>  inline double mysqrt<double>(double x){ return sqrt(x);}

template<typename X> inline X myabs(X x){ abort(); return sqrt(x);}
template<>  inline float myabs<float>(float x){ return fabsf(x);}
template<>  inline double myabs<double>(double x){ return fabs(x);}

//che cos'e' questo odore? - non era un sospiro.
inline Real heaviside(const Real phi, const Real inv_h) //11 FLOP
{
	const Real x = min((Real)1, max((Real)-1, phi * inv_h));
	const Real val_xneg = (((Real)-0.5)*x - 1)*x + ((Real)0.5);
	const Real val_xpos = (((Real)+0.5)*x - 1)*x + ((Real)0.5);
	const Real hs = x < 0 ? val_xneg : val_xpos;

	return hs;
}

inline Real reconstruct(const Real y0, const Real y1, const Real phi, const Real inv_h)//15 FLOP
{
	const Real hs = heaviside(phi, inv_h);

	return y0*hs + y1*(1-hs);
}

inline Real getgamma(const Real phi, const Real smoothlength, const Real gamma0, const Real gamma1)
{
	return reconstruct(gamma0, gamma1, phi, 1/smoothlength);
}

inline Real getPC(const Real phi, const Real smoothlength, const Real pc1, const Real pc2)
{
	return reconstruct(pc1, pc2, phi, 1/smoothlength);
}

#if defined(_QPX_) || defined(_QPXEMU_)

//QPX related functions
#ifdef _QPX_
#define _DIEGO_TRANSPOSE4(a, b, c, d)\
{\
const vector4double v01L = vec_perm(a, b, vec_gpci(00415));\
const vector4double v01H = vec_perm(a, b, vec_gpci(02637));\
const vector4double v23L = vec_perm(c, d, vec_gpci(00415));\
const vector4double v23H = vec_perm(c, d, vec_gpci(02637));\
\
a = vec_perm(v01L, v23L, vec_gpci(00145));\
b = vec_perm(v01L, v23L, vec_gpci(02367));\
c = vec_perm(v01H, v23H, vec_gpci(00145));\
d = vec_perm(v01H, v23H, vec_gpci(02367));\
}
#endif

template<int preclevel>
vector4double inline myreciprocal(vector4double a) {abort();}

#if defined(_BGQ_)
template<>
vector4double inline myreciprocal<1>(vector4double a)
{
  return vec_swdiv_nochk((const vector4double)(1), a);//vec_swdiv((vector4double)(1), a);
}
#endif

template<>
vector4double inline myreciprocal<0>(vector4double a)
{
	const vector4double Ra0 = vec_re(a);
	return vec_madd(vec_nmsub(Ra0, a, vec_splats(1)), Ra0, Ra0);
	//return vec_mul(Ra0, vec_nmsub(Ra0, a, vec_splats(2)));
	//return vec_nmsub(vec_mul(Ra0, a), Ra0, vec_add(Ra0, Ra0));
}

template<>
vector4double inline myreciprocal<2>(vector4double a)
{
  const vector4double Ra0 = vec_re(a);
  const vector4double Ra1= vec_madd(vec_nmsub(Ra0, a, vec_splats(1)), Ra0, Ra0);
  return vec_madd(vec_nmsub(Ra1, a, vec_splats(1)), Ra1, Ra1);
}

template<>
vector4double inline myreciprocal<-1>(vector4double a) { return vec_re(a); }

template< int preclevel>
int myreciprocal_flops()
{
	if (preclevel == 0) return 4;
	if (preclevel == -1) return 1;
	return 1; //we don't know if it is just 1
}

/*
template< int preclevel >
vector4double inline mydivision(vector4double a, vector4double b)
{
	return vec_swdiv(a, b);
}
*/
template<int preclevel>
vector4double inline mydivision(vector4double a, vector4double b)
{
	return vec_mul(a, myreciprocal<preclevel>(b));
}
/*
template<>
vector4double inline mydivision<-1>(vector4double a, vector4double b)
{
	return vec_mul(a, vec_re(b));
}
*/
template< int preclevel>
int mydivision_flops()
{
	if (preclevel <= 0) return myreciprocal_flops<preclevel>() + 1;
	return 1; //we don't know if it is just 1
}

template<int preclevel>
inline vector4double mysqrt(const vector4double a) {abort(); }

template<>
inline vector4double mysqrt<0>(const vector4double a)
{
	const vector4double invz =  vec_rsqrte(a);
	const vector4double z = vec_re(invz);
	const vector4double tmp = vec_msub(z, z, a);
	//const vector4double candidate = vec_nmsub(tmp, vec_mul(invz, vec_splats(0.5f)), z);
	const vector4double candidate = vec_nmsub(tmp, vec_madd(invz, vec_splats(-0.5f), invz), z);
	return candidate;//vec_sel(candidate, vec_splats(0.f), vec_neg(a));//is it necessary?
}

#if defined(_BGQ_)
template<>
inline vector4double mysqrt<1>(const vector4double a)
{
  return vec_swsqrt_nochk(a);
}
#endif

template<>
inline vector4double mysqrt<-1>(const vector4double a)
{
  return vec_mul(a, vec_rsqrte(a));
}

template<>
inline vector4double mysqrt<2>(const vector4double a)
{
const vector4double z = vec_splats(0);
const vector4double estimate = vec_rsqrte( a );
const vector4double estimateSquared = vec_madd( estimate, estimate, z );
const vector4double halfEstimate = vec_madd( estimate, vec_splats(0.5), z );
const vector4double candidate  = vec_madd(a,vec_madd(vec_nmsub( a, estimateSquared, vec_splats(1.0)), halfEstimate, estimate ), z);

 return candidate;//vec_sel(candidate, vec_splats(0.f), vec_neg(a));//is it necessary?
}

template< int preclevel>
int mysqrt_flops()
{
	return 9;
}

inline vector4double mymin(const vector4double a, const vector4double b)
{
#ifdef _QPXEMU_
	return _mm_min_ps(a,b);
#else
	return vec_sel(a, b, vec_cmpgt(a, b));
#endif
}

inline vector4double mymax(const vector4double a, const vector4double b)
{
#ifdef _QPXEMU_
	return _mm_max_ps(a,b);
#else
	return vec_sel(a, b, vec_cmplt(a, b));
#endif
}

inline int myminmax_flops()
{
	return 2;
}
#else

template< int preclevel>
int myreciprocal_flops()
{
	if (preclevel == 0) return 4;
	if (preclevel == -1) return 1;
	return 1; //we don't know if it is just 1
}

template< int preclevel>
int mysqrt_flops()
{
	return 9;
}

template< int preclevel>
int mydivision_flops()
{
	if (preclevel <= 0) return myreciprocal_flops<preclevel>() + 1;
	return 1; //we don't know if it is just 1
}

inline int myminmax_flops()
{
	return 2;
}
#endif


//SSE-related functions
#ifdef _QPXEMU_
#include <xmmintrin.h>
#ifdef __INTEL_COMPILER
inline __m128 operator+(__m128 a, __m128 b){ return _mm_add_ps(a, b); }
inline __m128 operator&(__m128 a, __m128 b){ return _mm_and_ps(a, b); }
inline __m128 operator|(__m128 a, __m128 b){ return _mm_or_ps(a, b); }
inline __m128 operator*(__m128 a, __m128 b){ return _mm_mul_ps(a, b); }
inline __m128 operator-(__m128 a,  __m128 b){ return _mm_sub_ps(a, b); }
inline __m128 operator/(__m128 a, __m128 b){ return _mm_div_ps(a, b); }
inline __m128d operator+(__m128d a, __m128d b){ return _mm_add_pd(a, b); }
inline __m128d operator*(__m128d a, __m128d b){ return _mm_mul_pd(a, b); }
inline __m128d operator-(__m128d a, __m128d b){ return _mm_sub_pd(a, b); }
inline __m128d operator/(__m128d a, __m128d b){ return _mm_div_pd(a, b); }
inline __m128d operator&(__m128d a, __m128d b){ return _mm_and_pd(a, b); }
inline __m128d operator|(__m128d a, __m128d b){ return _mm_or_pd(a, b); }
#endif
#endif

inline void printKernelName(string kernelname)
{
	cout << endl;
	for (int i=0; i<100; i++)
		cout << "=";
	cout << endl << endl;
	cout << "KERNEL - " << kernelname.c_str() << endl << endl;
}

inline void printEndKernelTest()
{
	cout << endl;
	for (int i=0; i<100; i++)
		cout << "=";
	cout << endl << endl;
}

inline void printAccuracyTitle()
{
	string acc = " ACCURACY ";
	int l = (80-acc.size())/2;
	cout << "\t";
	for (int i=0; i<l; i++)
		cout << "-";
	cout << acc.c_str();
	cout << "";
	for (int i=0; i<l; i++)
		cout << "-";
	cout << endl;
}

inline void printPerformanceTitle()
{
	string perf = " PERFORMANCE ";
	int l = (80-perf.size())/2;
	cout << "\t";
	for (int i=0; i<l; i++)
		cout << "-";
	cout << perf.c_str();
	cout << "";
	for (int i=0; i<80-l-(int)perf.size(); i++)
		cout << "-";
	cout << endl;
}

inline void printEndLine()
{
	cout << "\t";
	for (int i=0; i<80; i++)
		cout << "-";
	cout << endl;
}

inline void awkAcc(string kernelname,
				   double linf0, double linf1, double linf2, double linf3, double linf4, double linf5,
				   double l10, double l11, double l12, double l13, double l14, double l15)
{
	cout << "awkAccLinf\t" <<  kernelname.c_str();
	cout << "\t" << setprecision(4) << linf0;
	cout << "\t" << setprecision(4) << linf1;
	cout << "\t" << setprecision(4) << linf2;
	cout << "\t" << setprecision(4) << linf3;
	cout << "\t" << setprecision(4) << linf4;
	cout << "\t" << setprecision(4) << linf5;
	cout << endl;

	cout << "awkAccL1\t" <<  kernelname.c_str();
	cout << "\t" << setprecision(4) << l10;
	cout << "\t" << setprecision(4) << l11;
	cout << "\t" << setprecision(4) << l12;
	cout << "\t" << setprecision(4) << l13;
	cout << "\t" << setprecision(4) << l14;
	cout << "\t" << setprecision(4) << l15;
	cout << endl;
}
