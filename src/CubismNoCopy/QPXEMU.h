/*
 *  QPXEMU.h
 *
 *
 *  Created by Diego Rossinelli on 4/10/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#ifdef _QPXEMU_
#include <xmmintrin.h>
#define vector4double __m128
#define _DIEGO_TRANSPOSE4(a,b,c,d) _MM_TRANSPOSE4_PS(a,b,c,d)
#define vec_lda(a, b) _mm_load_ps(b + ((a)/4))
#define vec_ld(a, b) _mm_loadu_ps(b + ((a)/4))
#define vec_sta(a, b, c) _mm_store_ps(c + ((b)/4), a)
#define vec_st(a, b, c) _mm_storeu_ps(c + ((b)/4), a)
#define vec_mul(a, b) (a) * (b)
#define vec_sub(a, b) (a) - (b)
#define vec_add(a, b) (a) + (b)
#define vec_madd(a, b, c) (a) * (b) + (c)
#define vec_msub(a, b, c) (a) * (b) - (c)
#define vec_nmsub(a, b, c) (c) - (a) * (b)
#define vec_splats(a) _mm_set1_ps(a)
#define vec_swdiv(a, b) _mm_div_ps(a, b)
#define vec_res(a) _mm_rcp_ps(a)
#define vec_re(a) _mm_rcp_ps(a)
#define vec_rsqrtes(a) _mm_rsqrt_ps(a)
#define vec_rsqrte(a) _mm_rsqrt_ps(a)
#define vec_neg(a) _mm_setzero_ps() - (a)
#define vec_gpci(a) a
#define vec_perm(a,b,code) myshuffle<code>(a,b)
#define vec_extract(a, b) my_extract<b>(a)
#define vec_cmpgt(a, b) my_compare(a,b)
#define vec_cmpeq(a, b) my_cmpequal(a,b)
#define vec_min(a, b) _mm_min_ps(a, b)
#define vec_max(a, b) _mm_max_ps(a, b)

static const unsigned int absvalmask1 = 0x7fffffff;
const float * const absmaskptr = (float *)&absvalmask1;
const __m128 absmask = _mm_set1_ps(*absmaskptr);

inline __m128 vec_abs(__m128 a)
{
	return _mm_and_ps(absmask, a);
}

template<int index>
inline float my_extract(__m128 a)
{
	float x[4];

	_mm_storeu_ps(x, a);

	return x[index];
}

template<int>
inline __m128 myshuffle(__m128 a, __m128 b);

template<>
inline __m128 myshuffle<01114>(__m128 a, __m128 b)
{
	return _mm_shuffle_ps(a, _mm_shuffle_ps(a, b, _MM_SHUFFLE(0,0,1,1)), _MM_SHUFFLE(3,1,1,1));
}

template<>
inline __m128 myshuffle<2323>(__m128 a, __m128 b)
{
	return _mm_shuffle_ps(a, a, _MM_SHUFFLE(3,2,3,2));
}

template<>
inline __m128 myshuffle<01234>(__m128 a, __m128 b)
{
	return _mm_shuffle_ps(a, _mm_shuffle_ps(a,b, _MM_SHUFFLE(0,0,3,3)), _MM_SHUFFLE(3,0,2,1));
}

template<>
inline __m128 myshuffle<02345>(__m128 a, __m128 b)
{
	return _mm_shuffle_ps(a, b, _MM_SHUFFLE(1,0,3,2));
}

template<>
inline __m128 myshuffle<03456>(__m128 a, __m128 b)
{
	return _mm_shuffle_ps(_mm_shuffle_ps(a,b, _MM_SHUFFLE(0,0,3,3)), b, _MM_SHUFFLE(2,1,3,0));
}

inline __m128 vec_sel(const __m128 a, const __m128 b, const __m128 c)
{
	const __m128 bwins = _mm_cmpge_ps(c, _mm_setzero_ps());
	return _mm_or_ps(_mm_and_ps(bwins, b), _mm_andnot_ps(bwins, a));
}

inline __m128 my_compare(__m128 a, __m128 b)
{
    const __m128 res = _mm_cmpgt_ps(a,b);
	return  _mm_or_ps(_mm_and_ps(res,_mm_set1_ps(1.0)),_mm_andnot_ps(res,_mm_set1_ps(-1.0)));
}

inline __m128 my_cmpequal(__m128 a, __m128 b)
{
    const __m128 res = _mm_cmpeq_ps(a,b);
	return  _mm_or_ps(_mm_and_ps(res,_mm_set1_ps(1.0)),_mm_andnot_ps(res,_mm_set1_ps(-1.0)));
}

#define __align(_ALIGNBYTES_) __attribute__((__aligned__(_ALIGNBYTES_)))

#endif

