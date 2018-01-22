/*
 *  Diffusion_QPX.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 3/2/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */
#pragma once

//#include <xmmintrin.h>

#include "Diffusion_CPP.h"
#include "DivTensor_QPX.h"

class Diffusion_QPX: public virtual Diffusion_CPP, public virtual DivTensor_QPX
{
public:

 Diffusion_QPX(const Real a=1, const Real mu1 = 1, const Real mu2 = 2, const Real G1 = 1/(2.1-1), const Real G2 = 1/(2.1-1), const Real h = 1, const Real smooth_length=1, const Real dtinvh = 1):
	DivTensor_CPP(a, dtinvh, h, 1.0), //this is the base class of Diffusion_CPP and DivTensor_QPX
	  Diffusion_CPP(a, mu1, mu2, G1, G2, h, smooth_length, dtinvh), DivTensor_QPX(a, dtinvh, h, 1.0)
	{
	}

protected:

	inline vector4double better_rcp(const vector4double a)
	{
		const vector4double Ra0 = myreciprocal<preclevel>(a);
		return vec_sub(vec_add(Ra0, Ra0), vec_mul(vec_mul(Ra0, a), Ra0));
	}

	//used only for biphase
	inline vector4double _compute_mu(const vector4double G,  const vector4double G2,
							  const vector4double ls_factor, const vector4double mu1, const vector4double mu2, const vector4double one)
	{
		const vector4double x = (G-G2)*ls_factor;
		const vector4double lambda = vec_min(vec_max(x ,vec_splats(0.0)), one);

		return mu2*(one-lambda) + mu1*lambda;
	}

	template<bool biphase>
	inline void _convert(const float * const pt, const float ls_factor, float&u, float&v, float&w, float& mu)
	{
		if (biphase)
		{
			const float inv_rho = (float)(1.f/pt[0]);

			u = pt[1]*inv_rho;
			v = pt[2]*inv_rho;
			w = pt[3]*inv_rho;

			const float x = (pt[5]-G2)*ls_factor;
			const float lambda = min( max(x ,(float)0), (float)1);

			mu = mu2*((float)1-lambda) + mu1*lambda;
		}
		else
		{
			const float inv_rho = (float)(1.f/pt[0]);

			u = pt[1]*inv_rho;
			v = pt[2]*inv_rho;
			w = pt[3]*inv_rho;
			mu = mu1;
		}
	}

	template<bool biphase>
	void _convert_sse(const float * const gptfirst, const int gptfloats, const int rowgpts,
					  InputSOAf_ST& _u, InputSOAf_ST& _v, InputSOAf_ST& _w, InputSOAf_ST& _mu)
	{
		const float ls_factor = ((float)1/(G1 - G2));

		const vector4double F_1 = vec_splats(1);
		const vector4double F_G2 = vec_splats(G2);
		const vector4double F_LS = vec_splats(ls_factor);
		const vector4double F_MU1 = vec_splats(mu1);
		const vector4double F_MU2 = vec_splats(mu2);

		const int stride = gptfloats*rowgpts;
		static const int PITCHOUT = _u.PITCH;

		float * const ubase = &_u.ref(-1,-1);
		float * const vbase = &_v.ref(-1,-1);
		float * const wbase = &_w.ref(-1,-1);
		float * const mubase = &_mu.ref(-1,-1);

		for(int iy=0; iy<_BLOCKSIZE_+2; iy++)
		{
			float * const u = ubase + iy*PITCHOUT;
			float * const v = vbase + iy*PITCHOUT;
			float * const w = wbase + iy*PITCHOUT;
			float * const mu = mubase + iy*PITCHOUT;

			const float * const leftghost = gptfirst + stride*iy;

			_convert<biphase>(leftghost, ls_factor, u[0], v[0], w[0], mu[0]);

			for(int ix=1; ix<_BLOCKSIZE_+1; ix+=4)
			{
				const float * const gpt0 = leftghost + gptfloats*(ix+0);
				const float * const gpt1 = leftghost + gptfloats*(ix+1);
				const float * const gpt2 = leftghost + gptfloats*(ix+2);
				const float * const gpt3 = leftghost + gptfloats*(ix+3);

				/*static const int SL = 4;
				if (SL + iy <= _BLOCKSIZE_+1)
				{
				  _mm_prefetch((char *)(gpt0 + SL*rowgpts*gptfloats), _MM_HINT_T0);
				  _mm_prefetch((char *)(gpt1 + SL*rowgpts*gptfloats), _MM_HINT_T0);
				  _mm_prefetch((char *)(gpt2 + SL*rowgpts*gptfloats), _MM_HINT_T0);
				  _mm_prefetch((char *)(gpt3 + SL*rowgpts*gptfloats), _MM_HINT_T0);
				}*/

				vector4double data0 = vec_ld(0L,gpt0);
				vector4double data1 = vec_ld(0L,gpt1);
				vector4double data2 = vec_ld(0L,gpt2);
				vector4double data3 = vec_ld(0L,gpt3);

				_DIEGO_TRANSPOSE4(data0, data1, data2, data3);

#ifdef _PREC_DIV_
				const vector4double inv_rho = F_1/data0;
#else
				const vector4double inv_rho = better_rcp(data0);
#endif
				vec_sta(data1*inv_rho, 0L, u + ix);
				vec_sta(data2*inv_rho, 0L, v + ix);
				vec_sta(data3*inv_rho, 0L, w + ix);

				if (biphase)
                {
#ifdef _QPXEMU_
                    const __m128 dummy = _mm_set_ps(gpt3[5], gpt2[5], gpt1[5], gpt0[5]);
#else
                    const vector4double dummy = (vector4double)(gpt3[5], gpt2[5], gpt1[5], gpt0[5]);
#endif
					vec_sta(data0*_compute_mu(dummy, F_G2, F_LS, F_MU1, F_MU2, F_1), 0L, mu + ix);
                }
				else
					vec_sta(data0*F_MU1, 0L, mu + ix);
			}

			_convert<biphase>(leftghost + gptfloats*(_BLOCKSIZE_+1), ls_factor, u[_BLOCKSIZE_+1], v[_BLOCKSIZE_+1], w[_BLOCKSIZE_+1], mu[_BLOCKSIZE_+1]);
		}
	}

	void _convert(const float * const gptfirst, const int gptfloats, const int rowgpts)
	{
#ifndef _SP_COMP_
		printf("Diffusion_QPX::_convert: you should not be here in double precision. Aborting.\n");
		abort();
#else
		InputSOA_ST& u=ringu.ref(), &v=ringv.ref(), &w=ringw.ref(), &mu = ringls.ref();

		if (G1 == G2)
			_convert_sse<false>(gptfirst, gptfloats, rowgpts, u, v, w, mu);
		else
			_convert_sse<true>(gptfirst, gptfloats, rowgpts, u, v, w, mu);
#endif
}

#define	LEFT(W,C) \
myshuffle<03456>(W,C)

	void _xface(const InputSOAf_ST& _mu,
				const TempSOAf_ST& ux0, const TempSOAf_ST& uy0, const TempSOAf_ST& uz0, const TempSOAf_ST& ux1, const TempSOAf_ST& uy1, const TempSOA_ST& uz1,
				const TempSOAf_ST& vx0, const TempSOAf_ST& vy0, const TempSOAf_ST& vx1, const TempSOAf_ST& vy1,
				const TempSOAf_ST& wx0, const TempSOAf_ST& wz0, const TempSOAf_ST& wx1, const TempSOAf_ST& wz1)
	{
#ifndef _SP_COMP_
		printf("Diffusion_QPX::_xface: you should not be here in double precision. Aborting.\n");
		abort();
#else
		const float factor = (float)(2./3);

		_average_xface<false>(2-factor, ux0, ux1, txx);
		_average_xface<true>(-factor, vy0, vy1, txx);
		_average_xface<true>(-factor, wz0, wz1, txx);

		_average_xface<false>(1, vx0, vx1, txy);
		_average_xface<true>(1, uy0, uy1, txy);

		_average_xface<false>(1, uz0, uz1, txz);
		_average_xface<true>(1, wx0, wx1, txz);

		static const int PITCHIN = InputSOA_ST::PITCH;
		static const int PITCHOUT = TempPiXSOA_ST::PITCH;

		const float * const mubase = _mu.ptr(0,0);
		float * const txxbase = &txx.ref(0,0);
		float * const txybase = &txy.ref(0,0);
		float * const txzbase = &txz.ref(0,0);

		for(int iy=0; iy<TempPiXSOA_ST::NY; iy++)
		{
			const float * const muptr = mubase + iy*PITCHIN;
			float * const txxptr = txxbase + iy*PITCHOUT;
			float * const txyptr = txybase + iy*PITCHOUT;
			float * const txzptr = txzbase + iy*PITCHOUT;

			vector4double W =  vec_lda(0L,muptr - 4);

			for(int ix=0; ix<TempPiXSOA_ST::NX; ix+=4)
			{
				const vector4double C = vec_lda(0L,muptr + ix);
				const vector4double mu =  vec_splats(0.5)*(C + LEFT(W,C));

				vec_sta(mu*vec_lda(0L,txxptr + ix), 0L, txxptr + ix);
				vec_sta(mu*vec_lda(0L,txyptr + ix), 0L, txyptr + ix);
				vec_sta(mu*vec_lda(0L,txzptr + ix), 0L, txzptr + ix);

				W = C;
			}
		}
#endif
	}

	void _yface(const InputSOAf_ST& _mu,
				const TempSOAf_ST& ux0, const TempSOAf_ST& uy0, const TempSOAf_ST& ux1, const TempSOAf_ST& uy1,
				const TempSOAf_ST& vx0, const TempSOAf_ST& vy0, const TempSOAf_ST& vz0, const TempSOAf_ST& vx1, const TempSOAf_ST& vy1, const TempSOA_ST& vz1,
				const TempSOAf_ST& wy0, const TempSOAf_ST& wz0, const TempSOAf_ST& wy1, const TempSOAf_ST& wz1)
	{
#ifndef _SP_COMP_
		printf("Diffusion_QPX::_yface: you should not be here in double precision. Aborting.\n");
		abort();
#else
		const float factor = (float)(2./3);

		_average_yface<false>(1, uy0, uy1, tyx);
		_average_yface<true>(1, vx0, vx1, tyx);

		_average_yface<false>(-factor, ux0, ux1, tyy);
		_average_yface<true>(2-factor, vy0, vy1, tyy);
		_average_yface<true>(-factor, wz0, wz1, tyy);

		_average_yface<false>(1, vz0, vz1, tyz);
		_average_yface<true>(1, wy0, wy1, tyz);

		static const int PITCHIN = InputSOA_ST::PITCH;
		static const int PITCHOUT = TempPiYSOA_ST::PITCH;

		const float * const mubase = _mu.ptr(0,0);
		float * const tyxbase = &tyx.ref(0,0);
		float * const tyybase = &tyy.ref(0,0);
		float * const tyzbase = &tyz.ref(0,0);

		for(int iy=0; iy<TempPiYSOA_ST::NY; iy++)
		{
			const float * const muptr = mubase + iy*PITCHIN;
			float * const tyxptr = tyxbase + iy*PITCHOUT;
			float * const tyyptr = tyybase + iy*PITCHOUT;
			float * const tyzptr = tyzbase + iy*PITCHOUT;

			for(int ix=0; ix<TempPiYSOA_ST::NX; ix+=4)
			{
			  const vector4double mu =  vec_splats(0.5)*(vec_lda(0L,muptr + ix - PITCHIN) + vec_lda(0L,muptr + ix));

				vec_sta(mu*vec_lda(0L,tyxptr + ix), 0L, tyxptr + ix);
				vec_sta(mu*vec_lda(0L,tyyptr + ix), 0L, tyyptr + ix);
				vec_sta(mu*vec_lda(0L,tyzptr + ix), 0L, tyzptr + ix);
			}
		}
#endif
	}

	void _zface(const InputSOAf_ST& mu0, const InputSOAf_ST& mu1,
							   const TempSOAf_ST& ux, const TempSOAf_ST& uz,
							   const TempSOAf_ST& vy, const TempSOAf_ST& vz,
							   const TempSOAf_ST& wx, const TempSOAf_ST& wy, const TempSOAf_ST& wz,
							   TempPiZSOAf_ST& tzx, TempPiZSOAf_ST& tzy, TempPiZSOAf_ST& tzz)
	{
#ifndef _SP_COMP_
		printf("Diffusion_QPX::_zface: you should not be here in double precision. Aborting.\n");
		abort();
#else
		const float factor = (float)(2./3);

		_average_zface<false>(1, uz, tzx);
		_average_zface<true>(1, wx, tzx);

		_average_zface<false>(1, vz, tzy);
		_average_zface<true>(1, wy, tzy);

		_average_zface<false>(-factor, ux, tzz);
		_average_zface<true>(-factor, vy, tzz);
		_average_zface<true>(2-factor, wz, tzz);

		static const int PITCHIN = InputSOA_ST::PITCH;
		static const int PITCHOUT = TempPiZSOA_ST::PITCH;

		const float * const mubase0 = mu0.ptr(0,0);
		const float * const mubase1 = mu1.ptr(0,0);
		float * const tzxbase = &tzx.ref(0,0);
		float * const tzybase = &tzy.ref(0,0);
		float * const tzzbase = &tzz.ref(0,0);

		for(int iy=0; iy<TempPiZSOA_ST::NY; iy++)
		{
			const float * const muptr0 = mubase0 + iy*PITCHIN;
			const float * const muptr1 = mubase1 + iy*PITCHIN;

			float * const tzxptr = tzxbase + iy*PITCHOUT;
			float * const tzyptr = tzybase + iy*PITCHOUT;
			float * const tzzptr = tzzbase + iy*PITCHOUT;

			for(int ix=0; ix<TempPiZSOA_ST::NX; ix+=4)
			{
			  const vector4double mu = vec_splats(0.5)*(vec_lda(0L,muptr0 + ix) + vec_lda(0L,muptr1 + ix));

				vec_sta(mu*vec_lda(0L,tzxptr + ix), 0L, tzxptr + ix);
				vec_sta(mu*vec_lda(0L,tzyptr + ix), 0L, tzyptr + ix);
				vec_sta(mu*vec_lda(0L,tzzptr + ix), 0L, tzzptr + ix);
			}
		}
#endif
}

#ifndef _ALTERNATIVEIMPL_
	void _xmul(const InputSOAf_ST& _mu)
	{
#ifndef _SP_COMP_
		printf("Diffusion_QPX::_xmul: you should not be here in double precision. Aborting.\n");
		abort();
#else
		static const int PITCHIN = InputSOA_ST::PITCH;
		static const int PITCHOUT = TempPiXSOA_ST::PITCH;

		const float * const mubase = _mu.ptr(0,0);
		float * const txxbase = &txx.ref(0,0);
		float * const txybase = &txy.ref(0,0);
		float * const txzbase = &txz.ref(0,0);

		for(int iy=0; iy<TempPiXSOA_ST::NY; iy++)
		{
			const float * const muptr = mubase + iy*PITCHIN;
			float * const txxptr = txxbase + iy*PITCHOUT;
			float * const txyptr = txybase + iy*PITCHOUT;
			float * const txzptr = txzbase + iy*PITCHOUT;

			vector4double W =  vec_lda(0L,muptr - 4);

			for(int ix=0; ix<TempPiXSOA_ST::NX; ix+=4)
			{
				const vector4double C = vec_lda(0L,muptr + ix);
				const vector4double mu = vec_splats(0.5)*(C + LEFT(W,C));

				vec_sta(mu*vec_lda(0L,txxptr + ix), 0L, txxptr + ix);
				vec_sta(mu*vec_lda(0L,txyptr + ix), 0L, txyptr + ix);
				vec_sta(mu*vec_lda(0L,txzptr + ix), 0L, txzptr + ix);

				W = C;
			}
		}
#endif
	}

	void _ymul(const InputSOAf_ST& _mu)
	{
#ifndef _SP_COMP_
		printf("Diffusion_QPX::_ymul: you should not be here in double precision. Aborting.\n");
		abort();
#else

		static const int PITCHIN = InputSOA_ST::PITCH;
		static const int PITCHOUT = TempPiYSOA_ST::PITCH;

		const float * const mubase = _mu.ptr(0,0);
		float * const tyxbase = &tyx.ref(0,0);
		float * const tyybase = &tyy.ref(0,0);
		float * const tyzbase = &tyz.ref(0,0);

		for(int iy=0; iy<TempPiYSOA_ST::NY; iy++)
		{
			const float * const muptr = mubase + iy*PITCHIN;
			float * const tyxptr = tyxbase + iy*PITCHOUT;
			float * const tyyptr = tyybase + iy*PITCHOUT;
			float * const tyzptr = tyzbase + iy*PITCHOUT;

			for(int ix=0; ix<TempPiYSOA_ST::NX; ix+=4)
			{
			  const vector4double mu = vec_splats(0.5)*(vec_lda(0L,muptr + ix - PITCHIN) + vec_lda(0L,muptr + ix));

				vec_sta(mu*vec_lda(0L,tyxptr + ix), 0L, tyxptr + ix);
				vec_sta(mu*vec_lda(0L,tyyptr + ix), 0L, tyyptr + ix);
				vec_sta(mu*vec_lda(0L,tyzptr + ix), 0L, tyzptr + ix);
			}
		}
#endif
	}

	void _zmul(const InputSOAf_ST& mu0, const InputSOAf_ST& mu1, TempPiZSOAf_ST& tzx, TempPiZSOAf_ST& tzy, TempPiZSOAf_ST& tzz)
	{
#ifndef _SP_COMP_
		printf("Diffusion_QPX::_zmul: you should not be here in double precision. Aborting.\n");
		abort();
#else
		static const int PITCHIN = InputSOA_ST::PITCH;
		static const int PITCHOUT = TempPiZSOA_ST::PITCH;

		const float * const mubase0 = mu0.ptr(0,0);
		const float * const mubase1 = mu1.ptr(0,0);
		float * const tzxbase = &tzx.ref(0,0);
		float * const tzybase = &tzy.ref(0,0);
		float * const tzzbase = &tzz.ref(0,0);

		for(int iy=0; iy<TempPiZSOA_ST::NY; iy++)
		{
			const float * const muptr0 = mubase0 + iy*PITCHIN;
			const float * const muptr1 = mubase1 + iy*PITCHIN;

			float * const tzxptr = tzxbase + iy*PITCHOUT;
			float * const tzyptr = tzybase + iy*PITCHOUT;
			float * const tzzptr = tzzbase + iy*PITCHOUT;

			for(int ix=0; ix<TempPiZSOA_ST::NX; ix+=4)
			{
			  const vector4double mu = vec_splats(0.5)*(vec_lda(0L,muptr0 + ix) + vec_lda(0L,muptr1 + ix));

				vec_sta(mu*vec_lda(0L,tzxptr + ix), 0L, tzxptr + ix);
				vec_sta(mu*vec_lda(0L,tzyptr + ix), 0L, tzyptr + ix);
				vec_sta(mu*vec_lda(0L,tzzptr + ix), 0L, tzzptr + ix);
			}
		}
#endif
	}

public:

void compute(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
                 Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts)
  {
		_convert(srcfirst, srcfloats, rowsrcs);
		_input_next();

		_convert(srcfirst + srcfloats*slicesrcs, srcfloats, rowsrcs);

		const float factor = (float)(2./3);

		_corners(ringu(-1), ringu(0), ringnx.ref(), ringny.ref(), ringnz.ref());

		{
			const TempSOAf_ST& ux1 = ringnx.ref();
			const TempSOAf_ST& uz1 = ringnz.ref();

			_average_zface<false>(1, uz1, ringtzx.ref());
			_average_zface<false>(-factor, ux1, ringtzz.ref());
		}

		_corners(ringv(-1), ringv(0), gradvx.ref(), gradvy.ref(), gradvz.ref());

		{
			const TempSOAf_ST& vy1 = gradvy.ref();
			const TempSOAf_ST& vz1 = gradvz.ref();

			_average_zface<false>(1, vz1, ringtzy.ref());
			_average_zface<true>(-factor, vy1, ringtzz.ref());
		}

		_corners(ringw(-1), ringw(0), gradwx.ref(), gradwy.ref(), gradwz.ref());
		{
			const TempSOAf_ST& wx1 = gradwx.ref();
			const TempSOAf_ST& wy1 = gradwy.ref();
			const TempSOAf_ST& wz1 = gradwz.ref();

			_average_zface<true>(1, wx1, ringtzx.ref());
			_average_zface<true>(1, wy1, ringtzy.ref());
			_average_zface<true>(2-factor, wz1, ringtzz.ref());
		}

		_zmul(ringls(-1), ringls(), ringtzx.ref(), ringtzy.ref(), ringtzz.ref());

		_udot_tz(ringu(-1), ringv(-1), ringw(-1),  ringu(), ringv(), ringw(), ringtzx(), ringtzy(), ringtzz(), ringutz.ref());

		for(int islice=0; islice<_BLOCKSIZE_; islice++)
		{
			_tensors_next();
			_input_next();
			_grad_next();

			_convert(srcfirst + (islice+2)*srcfloats*slicesrcs, srcfloats, rowsrcs);

			_corners(ringu(-1), ringu(0), ringnx.ref(), ringny.ref(), ringnz.ref());

			{
				const TempSOAf_ST& ux0 = ringnx.ref(-1);
				const TempSOAf_ST& ux1 = ringnx.ref();
				const TempSOAf_ST& uy0 = ringny.ref(-1);
				const TempSOAf_ST& uy1 = ringny.ref();
				const TempSOAf_ST& uz0 = ringnz.ref(-1);
				const TempSOAf_ST& uz1 = ringnz.ref();

				_average_xface<false>(2-factor, ux0, ux1, txx);
				_average_yface<false>(-factor, ux0, ux1, tyy);
				_average_xface<false>(1, uy0, uy1, txy);
				_average_yface<false>(1, uy0, uy1, tyx);
				_average_xface<false>(1, uz0, uz1, txz);
				_average_zface<false>(1, uz1, ringtzx.ref());
				_average_zface<false>(-factor, ux1, ringtzz.ref());
			}

			_corners(ringv(-1), ringv(0), gradvx.ref(), gradvy.ref(), gradvz.ref());

			{
				const TempSOAf_ST& vx0 = gradvx.ref(-1);
				const TempSOAf_ST& vx1 = gradvx.ref();
				const TempSOAf_ST& vy0 = gradvy.ref(-1);
				const TempSOAf_ST& vy1 = gradvy.ref();
				const TempSOAf_ST& vz0 = gradvz.ref(-1);
				const TempSOAf_ST& vz1 = gradvz.ref();

				_average_xface<true>(-factor, vy0, vy1, txx);
				_average_xface<true>(1, vx0, vx1, txy);
				_average_yface<true>(1, vx0, vx1, tyx);
				_average_yface<true>(2-factor, vy0, vy1, tyy);
				_average_zface<true>(-factor, vy1, ringtzz.ref());
				_average_yface<false>(1, vz0, vz1, tyz);
				_average_zface<false>(1, vz1, ringtzy.ref());
			}

			_corners(ringw(-1), ringw(0), gradwx.ref(), gradwy.ref(), gradwz.ref());

			{
				const TempSOAf_ST& wx0 = gradwx.ref(-1);
				const TempSOAf_ST& wx1 = gradwx.ref();
				const TempSOAf_ST& wy0 = gradwy.ref(-1);
				const TempSOAf_ST& wy1 = gradwy.ref();
				const TempSOAf_ST& wz0 = gradwz.ref(-1);
				const TempSOAf_ST& wz1 = gradwz.ref();

				_average_xface<true>(1, wx0, wx1, txz);
				_average_zface<true>(1, wx1, ringtzx.ref());
				_average_yface<true>(1, wy0, wy1, tyz);
				_average_zface<true>(1, wy1, ringtzy.ref());
				_average_xface<true>(-factor, wz0, wz1, txx);
				_average_yface<true>(-factor, wz0, wz1, tyy);
				_average_zface<true>(2-factor, wz1, ringtzz.ref());
			}

			_xmul(ringls(-1));
			_ymul(ringls(-1));
			_zmul(ringls(-1), ringls(), ringtzx.ref(), ringtzy.ref(), ringtzz.ref());

			_udot_tx(ringu(-1), ringv(-1), ringw(-1));
			_udot_ty(ringu(-1), ringv(-1), ringw(-1));
			_udot_tz(ringu(-1), ringv(-1), ringw(-1), ringu(0), ringv(0), ringw(0), ringtzx(), ringtzy(), ringtzz(), ringutz.ref());

			_div_dxy();
			_div_dz(ringtzx(-1), ringtzy(-1), ringtzz(-1), ringutz(-1), ringtzx(0), ringtzy(0), ringtzz(0), ringutz(0));

			_copyback(dstfirst + islice*dstfloats*slicedsts, dstfloats, rowdsts);
		}
		}
#endif

#undef LEFT
};
