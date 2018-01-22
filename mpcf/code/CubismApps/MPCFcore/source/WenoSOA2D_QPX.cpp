/*
 **  WenoSOA2D_QPX.cpp
 **
 **
 **  Created by Diego Rossinelli on 2/12/13.
 **  Copyright 2013 ETH Zurich. All rights reserved.
 **
 **/

#include "common.h"
#include "WenoSOA2D_QPX.h"

#ifdef _WENO3_
#include "../../MPCFthread/source/Weno_QPX_3rdOrder.h"
#else
// #include "../../MPCFthread/source/WenoFused_QPX.h"
#include "../../MPCFthread/source/WenoFused_IS2_QPX.h"
#endif

#ifndef _MICROFUSION_
#define _MICROFUSION_ 2
#endif

//=============================== KERNELS ===========================================

struct __attribute__((__aligned__(_ALIGNBYTES_))) WenoScratchPad { Real tmp[TempSOA::NX][4];};

template<typename Tomato> inline void _qpx_xweno_minus(Real * const in, Real * const out)
{
	enum {
		NX = TempSOA::NX,
		NY = TempSOA::NY,
		INSTRIDE = InputSOA::PITCH,
		OUTSTRIDE = TempSOA::PITCH
	};

	Tomato weno_ftor;

	for(int dy=0; dy<NY; ++dy)
	{
		Real * const row = in + INSTRIDE * dy;

		/* optionally we could skip the reading of W4 and C4 */
#pragma unroll(1)
		for(int dx=0; dx<NX; dx += 4)
		{
			const vector4double W4 = vec_lda(0 * sizeof(Real), row + dx);
			const vector4double C4 = vec_lda(4 * sizeof(Real), row + dx);
			const vector4double E4 = vec_lda(8 * sizeof(Real), row + dx);

			const vector4double a = vec_perm(W4, C4, vec_gpci(01234));
			const vector4double b = vec_perm(W4, C4, vec_gpci(02345));
			const vector4double c = vec_perm(W4, C4, vec_gpci(03456));
			const vector4double d = C4;
			const vector4double e = vec_perm(C4, E4, vec_gpci(01234));

			const vector4double reconstruction = weno_ftor._weno(a, b, c, d, e);

			vec_sta(reconstruction, sizeof(Real) * (dx + OUTSTRIDE * dy), out);
		}
	}
}

template<typename Tomato> inline void _qpx_xweno_plus(Real * const in, Real * const out)
{
	enum {
		NX = TempSOA::NX,
		NY = TempSOA::NY,
		INSTRIDE = InputSOA::PITCH,
		OUTSTRIDE = TempSOA::PITCH
	};

	Tomato weno_ftor;

	for(int dy=0; dy<NY; ++dy)
	{

		Real * const row = in + INSTRIDE * dy;

		/* optionally we could skip the reading of W4 and C4 */
#pragma unroll(1)
		for(int dx=0; dx<NX; dx += 4)
		{
			const vector4double W4 = vec_lda(0 * sizeof(Real), row + dx);
			const vector4double C4 = vec_lda(4 * sizeof(Real), row + dx);
			const vector4double E4 = vec_lda(8 * sizeof(Real), row + dx);

			const vector4double a = vec_perm(W4, C4, vec_gpci(02345));
			const vector4double b = vec_perm(W4, C4, vec_gpci(03456));
			const vector4double c = C4;
			const vector4double d = vec_perm(C4, E4, vec_gpci(01234));
			const vector4double e = vec_perm(C4, E4, vec_gpci(02345));

			const vector4double reconstruction = weno_ftor._weno(a, b, c, d, e);

			vec_sta(reconstruction, sizeof(Real) * (dx + OUTSTRIDE * dy), out);
		}
	}
}

inline void _qpx_xweno_fused(Real * const in, Real * const outm, Real * const outp)
{
	enum {
		NX = TempSOA::NX,
		NY = TempSOA::NY,
		INSTRIDE = InputSOA::PITCH,
		OUTSTRIDE = TempSOA::PITCH
	};

#if _MICROFUSION_ == 1
	SoupMinus wminus;
	SoupPlus wplus;
#else // this is for  _MICROFUSION_ == 2
	Weno_QPX_fused fusedweno;
#endif
	for(int dy=0; dy<NY; ++dy)
	{
		Real * const row = in + INSTRIDE * dy;

#pragma unroll(1)
		for(int dx=0; dx<NX; dx += 4)
		{
			const vector4double W4 = vec_lda(0 * sizeof(Real), row + dx);
			const vector4double C4 = vec_lda(4 * sizeof(Real), row + dx);
			const vector4double E4 = vec_lda(8 * sizeof(Real), row + dx);

			const vector4double a = vec_perm(W4, C4, vec_gpci(01234));
			const vector4double b = vec_perm(W4, C4, vec_gpci(02345));
			const vector4double c = vec_perm(W4, C4, vec_gpci(03456));
			const vector4double d = C4;
			const vector4double e = vec_perm(C4, E4, vec_gpci(01234));
			const vector4double f = vec_perm(C4, E4, vec_gpci(02345));
#if _MICROFUSION_ == 1
			const vector4double reconstructionM = wminus._weno(a, b, c, d, e);
			const vector4double reconstructionP = wplus._weno(b, c, d, e, f);

			vec_sta(reconstructionM, sizeof(Real) * (dx + OUTSTRIDE * dy), outm);
			vec_sta(reconstructionP, sizeof(Real) * (dx + OUTSTRIDE * dy), outp);
#else // this is for  _MICROFUSION_ == 2
			fusedweno.weno_minus_plus_fused_opt2(a,b,c,d,e,f, outm + dx + OUTSTRIDE * dy, outp + dx + OUTSTRIDE * dy);
#endif
		}
	}
}


template<typename Chili, int start> inline void _qpx_yweno(Real * const in, Real * const out)
{
	enum {
		NX = TempSOA::NX,
		NY = TempSOA::NY,
		INSTRIDE = InputSOA::PITCH,
		OUTSTRIDE = TempSOA::PITCH
	};

	WenoScratchPad scratchpad;

	Real * const tmp = (Real *)&scratchpad.tmp[0][0];

	Chili weno_ftor;

	for(int dy=0; dy<NY; dy+=4)
	{
		Real * const ptr = (start * INSTRIDE) + &in[dy];

#pragma unroll(1)
		for(int dx=0; dx<NX; ++dx)
		{
			Real * const entry = ptr + dx * INSTRIDE;

			const vector4double a = vec_lda(0L, entry);
			const vector4double b = vec_lda(sizeof(Real) * INSTRIDE, entry);
			const vector4double c = vec_lda(sizeof(Real) * 2 * INSTRIDE, entry);
			const vector4double d = vec_lda(sizeof(Real) * 3 * INSTRIDE, entry);
			const vector4double e = vec_lda(sizeof(Real) * 4 * INSTRIDE, entry);

			const vector4double reconstruction = weno_ftor._weno(a, b, c, d, e);

			vec_sta(reconstruction, 0L, tmp + 4 * dx);
		}

#pragma unroll(1)
		for(int dx=0; dx<NX-1; dx += 4)
		{
			Real * const entry_in = tmp + 4 * dx;

			vector4double data0 = vec_lda(0L, entry_in);
			vector4double data1 = vec_lda(sizeof(Real) * 4, entry_in);
			vector4double data2 = vec_lda(sizeof(Real) * 4 * 2, entry_in);
			vector4double data3 = vec_lda(sizeof(Real) * 4 * 3, entry_in);

			_DIEGO_TRANSPOSE4(data0, data1, data2, data3);

			Real * const entry_out = out + dx + OUTSTRIDE * dy;

			vec_sta(data0, 0L, entry_out);
			vec_sta(data1, sizeof(Real) * OUTSTRIDE, entry_out);
			vec_sta(data2, sizeof(Real) * OUTSTRIDE * 2, entry_out);
			vec_sta(data3, sizeof(Real) * OUTSTRIDE * 3, entry_out);
		}

		{
			enum {
				A0 = NX-1,
				A1 = NX-1 + OUTSTRIDE,
				A2 = NX-1 + OUTSTRIDE * 2,
				A3 = NX-1 + OUTSTRIDE * 3,
				B = OUTSTRIDE
			};

			const int base = B * dy;

			out[A0 + base] = scratchpad.tmp[NX-1][0];
			out[A1 + base] = scratchpad.tmp[NX-1][1];
			out[A2 + base] = scratchpad.tmp[NX-1][2];
			out[A3 + base] = scratchpad.tmp[NX-1][3];
		}
	}
}

inline void _qpx_yweno_fused(Real * const in, Real * const outM, Real * const outP)
{
	enum {
		NX = TempSOA::NX,
		NY = TempSOA::NY,
		INSTRIDE = InputSOA::PITCH,
		OUTSTRIDE = TempSOA::PITCH
	};

	WenoScratchPad scratchpadM, scratchpadP;

	Real * const tmpM = (Real *)&scratchpadM.tmp[0][0];
	Real * const tmpP = (Real *)&scratchpadP.tmp[0][0];

#if _MICROFUSION_ == 1
	SoupMinus wminus;
	SoupPlus wplus;
#else // this is for  _MICROFUSION_ == 2
	Weno_QPX_fused fusedweno;
#endif
	for(int dy=0; dy<NY; dy+=4)
	{
		Real * const ptr = &in[dy];

#pragma unroll(1)
		for(int dx=0; dx<NX; ++dx)
		{
			Real * const entry = ptr + dx * INSTRIDE;

			const vector4double a = vec_lda(0L, entry);
			const vector4double b = vec_lda(sizeof(Real) * INSTRIDE, entry);
			const vector4double c = vec_lda(sizeof(Real) * 2 * INSTRIDE, entry);
			const vector4double d = vec_lda(sizeof(Real) * 3 * INSTRIDE, entry);
			const vector4double e = vec_lda(sizeof(Real) * 4 * INSTRIDE, entry);
			const vector4double f = vec_lda(sizeof(Real) * 5 * INSTRIDE, entry);
#if _MICROFUSION_ == 1
			const vector4double reconstructionM = wminus._weno(a, b, c, d, e);
			const vector4double reconstructionP = wplus._weno(b, c, d, e, f);

			vec_sta(reconstructionM, 0L, tmpM + 4 * dx);
			vec_sta(reconstructionP, 0L, tmpP + 4 * dx);
#else // this is for  _MICROFUSION_ == 2
			fusedweno.weno_minus_plus_fused_opt2(a,b,c,d,e,f, tmpM + 4 * dx, tmpP + 4 * dx);
#endif
		}


#pragma unroll(1)
		for(int dx=0; dx<NX-1; dx += 4)
		{
			{
				Real * const entry_in = tmpM + 4 * dx;

				vector4double data0 = vec_lda(0L, entry_in);
				vector4double data1 = vec_lda(sizeof(Real) * 4, entry_in);
				vector4double data2 = vec_lda(sizeof(Real) * 4 * 2, entry_in);
				vector4double data3 = vec_lda(sizeof(Real) * 4 * 3, entry_in);

				_DIEGO_TRANSPOSE4(data0, data1, data2, data3);

				Real * const entry_out = outM + dx + OUTSTRIDE * dy;

				vec_sta(data0, 0L, entry_out);
				vec_sta(data1, sizeof(Real) * OUTSTRIDE, entry_out);
				vec_sta(data2, sizeof(Real) * OUTSTRIDE * 2, entry_out);
				vec_sta(data3, sizeof(Real) * OUTSTRIDE * 3, entry_out);
			}

			{
				Real * const entry_in = tmpP + 4 * dx;

				vector4double data0 = vec_lda(0L, entry_in);
				vector4double data1 = vec_lda(sizeof(Real) * 4, entry_in);
				vector4double data2 = vec_lda(sizeof(Real) * 4 * 2, entry_in);
				vector4double data3 = vec_lda(sizeof(Real) * 4 * 3, entry_in);

				_DIEGO_TRANSPOSE4(data0, data1, data2, data3);

				Real * const entry_out = outP + dx + OUTSTRIDE * dy;

				vec_sta(data0, 0L, entry_out);
				vec_sta(data1, sizeof(Real) * OUTSTRIDE, entry_out);
				vec_sta(data2, sizeof(Real) * OUTSTRIDE * 2, entry_out);
				vec_sta(data3, sizeof(Real) * OUTSTRIDE * 3, entry_out);
			}
		}

		{
			enum {
				A0 = NX-1,
				A1 = NX-1 + OUTSTRIDE,
				A2 = NX-1 + OUTSTRIDE * 2,
				A3 = NX-1 + OUTSTRIDE * 3,
				B = OUTSTRIDE
			};

			const int base = B * dy;

			outM[A0 + base] = scratchpadM.tmp[NX-1][0];
			outM[A1 + base] = scratchpadM.tmp[NX-1][1];
			outM[A2 + base] = scratchpadM.tmp[NX-1][2];
			outM[A3 + base] = scratchpadM.tmp[NX-1][3];

			outP[A0 + base] = scratchpadP.tmp[NX-1][0];
			outP[A1 + base] = scratchpadP.tmp[NX-1][1];
			outP[A2 + base] = scratchpadP.tmp[NX-1][2];
			outP[A3 + base] = scratchpadP.tmp[NX-1][3];
		}
	}
}


template<typename Potato> inline void _qpx_zweno(Real * const a, Real * const b,
												 Real * const c, Real * const d,
												 Real * const e , Real * const out)
{
	enum {
		NX = TempSOA::NX,
		NY = TempSOA::NY,
		INSTRIDE = InputSOA::PITCH,
		OUTSTRIDE = TempSOA::PITCH
	};

	WenoQPX<Potato> wenoqpx;

	for(int dy=0; dy<NY; dy++)
		wenoqpx.stride<NX>(a + INSTRIDE * dy, b + INSTRIDE * dy, c + INSTRIDE * dy,
						   d + INSTRIDE * dy, e + INSTRIDE * dy, out + OUTSTRIDE * dy);
}

inline void _qpx_zweno_fused(Real * const _a, Real * const _b, Real * const _c, Real * const _d, Real * const _e, Real * const _f,
							 Real * const _outM, Real * const _outP)
{
	enum {
		NX = TempSOA::NX,
		NY = TempSOA::NY,
		INSTRIDE = InputSOA::PITCH,
		OUTSTRIDE = TempSOA::PITCH
	};

#if _MICROFUSION_ == 1
	SoupMinus wminus;
	SoupPlus wplus;
#else // this is for  _MICROFUSION_ == 2
	Weno_QPX_fused fusedweno;
#endif
	for(int dy=0; dy<NY; dy++)
	{
		Real * const a = _a + INSTRIDE * dy;
		Real * const b = _b + INSTRIDE * dy;
		Real * const c = _c + INSTRIDE * dy;
		Real * const d = _d + INSTRIDE * dy;
		Real * const e = _e + INSTRIDE * dy;
		Real * const f = _f + INSTRIDE * dy;

		Real * const outM = _outM + OUTSTRIDE * dy;
		Real * const outP = _outP + OUTSTRIDE * dy;

		for(int i=0; i<NX; i += 4)
		{
			const vector4double mya = vec_lda(0L, a + i);
			const vector4double myb = vec_lda(0L, b + i);
			const vector4double myc = vec_lda(0L, c + i);
			const vector4double myd = vec_lda(0L, d + i);
			const vector4double mye = vec_lda(0L, e + i);
			const vector4double myf = vec_lda(0L, f + i);

#if _MICROFUSION_ == 1
			const vector4double resultM = wminus._weno(mya, myb, myc, myd, mye);
			const vector4double resultP = wplus._weno(myb, myc, myd, mye, myf);

			vec_sta(resultM, 0L, outM + i);
			vec_sta(resultP, 0L, outP + i);
#else // this is for  _MICROFUSION_ == 2
			fusedweno.weno_minus_plus_fused_opt2(mya,myb,myc,myd,mye,myf, outM + i, outP + i);
#endif
		}
	}
}

//=============================== DRIVERS ===========================================

void WenoSOA2D_QPX::xcompute(const InputSOA& in, TempSOA& outm, TempSOA& outp) const
{
#if _MICROFUSION_ == 0
	_qpx_xweno_minus<SoupMinus>(const_cast<Real *>(in.ptr(-4,0)), &outm.ref(0,0));
	_qpx_xweno_plus<SoupPlus>(const_cast<Real *>(in.ptr(-4,0)), &outp.ref(0,0));
#else // this is for  _MICROFUSION_ == 1 OR 2
	_qpx_xweno_fused(const_cast<Real *>(in.ptr(-4,0)), &outm.ref(0,0), &outp.ref(0,0));
#endif
}

void WenoSOA2D_QPX::ycompute(const InputSOA& in, TempSOA& outm, TempSOA& outp) const

{
#if _MICROFUSION_ == 0
	_qpx_yweno<SoupMinus, 0>(const_cast<Real *>(in.ptr(0,-3)), &outm.ref(0,0));
	_qpx_yweno<SoupPlus, 1>(const_cast<Real *>(in.ptr(0,-3)), &outp.ref(0,0));
#else // this is for  _MICROFUSION_ == 1 OR 2
	_qpx_yweno_fused(const_cast<Real *>(in.ptr(0,-3)), &outm.ref(0,0), &outp.ref(0,0));
#endif
}

void WenoSOA2D_QPX::zcompute(const int r, const RingInputSOA& in, TempSOA& outm, TempSOA& outp) const

{
#if _MICROFUSION_ == 0
	_qpx_zweno<SoupMinus>(const_cast<Real *>(in(r-3).ptr(0,0)),
						  const_cast<Real *>(in(r-2).ptr(0,0)),
						  const_cast<Real *>(in(r-1).ptr(0,0)),
						  const_cast<Real *>(in(r).ptr(0,0)),
						  const_cast<Real *>(in(r+1).ptr(0,0)),
						  &outm.ref(0,0));

	_qpx_zweno<SoupPlus>(const_cast<Real *>(in(r-2).ptr(0,0)),
						 const_cast<Real *>(in(r-1).ptr(0,0)),
						 const_cast<Real *>(in(r).ptr(0,0)),
						 const_cast<Real *>(in(r+1).ptr(0,0)),
						 const_cast<Real *>(in(r+2).ptr(0,0)), &outp.ref(0,0));
#else // this is for  _MICROFUSION_ == 1 OR 2
	_qpx_zweno_fused(const_cast<Real *>(in(r-3).ptr(0,0)),
			const_cast<Real *>(in(r-2).ptr(0,0)),
			const_cast<Real *>(in(r-1).ptr(0,0)),
			const_cast<Real *>(in(r).ptr(0,0)),
			const_cast<Real *>(in(r+1).ptr(0,0)),
			const_cast<Real *>(in(r+2).ptr(0,0)), &outm.ref(0,0), &outp.ref(0,0));
#endif
}

