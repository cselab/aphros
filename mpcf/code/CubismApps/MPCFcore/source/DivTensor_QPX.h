/*
 *  DivTensor_QPX.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 3/2/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */
#pragma once

//#include <xmmintrin.h>
#include "DivTensor_CPP.h"

class DivTensor_QPX: public virtual DivTensor_CPP
{	
#define	LEFT(W,C) \
myshuffle<03456>(W,C)
	
#define RIGHT(C, E)	\
myshuffle<01234>(C,E)
	
protected:
	//virtual method of DivTensor_CPP overridden here
      void _corners(const InputSOAf_ST& _ls0, const InputSOAf_ST& _ls1, TempSOAf_ST& _nx, TempSOAf_ST& _ny, TempSOAf_ST& _nz)
	{	
		static const int PITCHIN = _ls0.PITCH;
		const float * const ls0base = _ls0.ptr(0,0);
		const float * const ls1base = _ls1.ptr(0,0);
		
		static const int PITCHOUT = _nx.PITCH;
		float * const nxbase = &_nx.ref(0,0);
		float * const nybase = &_ny.ref(0,0);
		float * const nzbase = &_nz.ref(0,0);
		
		for(int iy=0; iy<TempSOA_ST::NY; iy++)
		{
			const float * const ls0 = ls0base + iy*PITCHIN;
			const float * const ls1 = ls1base + iy*PITCHIN;
			
			float * const nx = nxbase + iy*PITCHOUT;
			float * const ny = nybase + iy*PITCHOUT;
			float * const nz = nzbase + iy*PITCHOUT;
			
			vector4double WA0 = vec_lda(0L,ls0 - 4);
			vector4double WB0 = vec_lda(0L,ls0 - PITCHIN -4);
			vector4double WA1 = vec_lda(0L,ls1 - 4);
			vector4double WB1 = vec_lda(0L,ls1 - PITCHIN -4);
			
			for(int ix=0; ix<TempSOA_ST::NX; ix+=4)
			{	
				const vector4double c00m1 = vec_lda(0L,ls0 + ix);
				const vector4double cm10m1 = LEFT(WA0, c00m1);
				const vector4double c0m1m1 = vec_lda(0L,ls0 + ix - PITCHIN);
				const vector4double cm1m1m1 = LEFT(WB0, c0m1m1);
				
				const vector4double c000 = vec_lda(0L,ls1 + ix);
				const vector4double cm100 = LEFT(WA1, c000);
				const vector4double c0m10 = vec_lda(0L,ls1 + ix - PITCHIN);
				const vector4double cm1m10 = LEFT(WB1, c0m10);
				
				vec_sta(c000-cm100 + c0m10-cm1m10 + c00m1-cm10m1 + c0m1m1-cm1m1m1, 0L, nx + ix);
				vec_sta(cm100-cm1m10 + c000-c0m10 + cm10m1-cm1m1m1  + c00m1-c0m1m1, 0L, ny + ix);
				vec_sta(cm1m10-cm1m1m1 + cm100-cm10m1 + c0m10-c0m1m1  + c000-c00m1, 0L, nz + ix);
				
				WA0 = c00m1;
				WB0 = c0m1m1;
				WA1 = c000;
				WB1 = c0m10;
			}
		}
	}
	
	void _copyback(float * const gptfirst, const int gptfloats, const int rowgpts)
	{	
#ifndef _SP_COMP_
		printf("DivTensor_QPX::_copyback: you should not be here in double precision. Aborting.\n");
		abort();
#else
		const vector4double F = vec_splats(sigma*dtinvh/(16*h));
		const vector4double A = vec_splats(a);
		const vector4double F2= F*vec_splats((float)0.5);
		
		static const int PITCHIN = OutputSOA::PITCH;
		const int stride = gptfloats*rowgpts;
		
		const float * const baserhsu = rhsu.ptr(0,0);
		const float * const baserhsv = rhsv.ptr(0,0);
		const float * const baserhsw = rhsw.ptr(0,0);
		const float * const baserhss = rhss.ptr(0,0);
		
		for(int iy=0; iy<OutputSOA::NY; iy++)
		{
			const float * const ptrrhsu = baserhsu + iy*PITCHIN;
			const float * const ptrrhsv = baserhsv + iy*PITCHIN;
			const float * const ptrrhsw = baserhsw + iy*PITCHIN;
			const float * const ptrrhss = baserhss + iy*PITCHIN;
			
			float * const basegp = 1 + gptfirst + iy*stride;
			
			for(int ix=0; ix<OutputSOA::NX; ix+=4)
			{
				float * const gp0 = basegp + gptfloats*(ix+0);
				float * const gp1 = basegp + gptfloats*(ix+1);
				float * const gp2 = basegp + gptfloats*(ix+2);
				float * const gp3 = basegp + gptfloats*(ix+3);
				
				vector4double rhs0 = F*vec_lda(0L,ptrrhsu + ix);
				vector4double rhs1 = F*vec_lda(0L,ptrrhsv + ix);
				vector4double rhs2 = F*vec_lda(0L,ptrrhsw + ix);
				vector4double rhs3 = F2*vec_lda(0L,ptrrhss + ix);
				
				_DIEGO_TRANSPOSE4(rhs0, rhs1, rhs2, rhs3);
				
				vec_st(A*vec_ld(0L,gp0) + rhs0, 0L, gp0);
				vec_st(A*vec_ld(0L,gp1) + rhs1, 0L, gp1);
				vec_st(A*vec_ld(0L,gp2) + rhs2, 0L, gp2);
				vec_st(A*vec_ld(0L,gp3) + rhs3, 0L, gp3);
			}
		}
#endif
	}
  
	template<int PITCHIN0, int PITCHIN1, int PITCHOUT>
	void _div_dxy_sse(const float * const basesrc0, const float * const basesrc1, float * const basedest)
	{
		for(int iy=0; iy<OutputSOA::NY; iy++)
		{
			const float * const src0 = basesrc0 + iy*PITCHIN0;
			const float * const src1 = basesrc1 + iy*PITCHIN1;
			float * const dest = basedest + iy*PITCHOUT;
			
			vector4double C = vec_lda(0L,src0);
			
			for(int ix=0; ix<OutputSOA::NX; ix+=4)
			{
				const vector4double E = vec_lda(0L,src0 + ix + 4);
				
				vec_sta(RIGHT(C, E) - C +
							 vec_lda(0L,src1 + ix + PITCHIN1) - vec_lda(0L,src1 + ix), 0L, dest + ix);
							 
				C = E;
			}
		}
	}

	template<int PITCHIN, int PITCHOUT>
	void _div_dz_sse(const float * const basesrc0, const float * const basesrc1, float * const basedest)
	{
		for(int iy=0; iy<OutputSOA::NY; iy++)
		{
			const float * const src0 = basesrc0 + iy*PITCHIN;
			const float * const src1 = basesrc1 + iy*PITCHIN;
			float * const dest = basedest + iy*PITCHOUT;
			
			for(int ix=0; ix<OutputSOA::NX; ix+=4)
				vec_sta(vec_lda(0L,dest + ix) + vec_lda(0L,src0 + ix) - vec_lda(0L,src1 + ix), 0L, dest + ix);
		}
	}
	
#ifdef _SP_COMP_
	void _div_dxy()
	{	
		static const int PIN0 = TempPiXSOA_ST::PITCH;
		static const int PIN1 = TempPiYSOA_ST::PITCH;
		static const int POUT = OutputSOA::PITCH;
		
		_div_dxy_sse<PIN0, PIN1, POUT>(txx.ptr(0,0), tyx.ptr(0,0), &rhsu.ref(0,0));
		_div_dxy_sse<PIN0, PIN1, POUT>(txy.ptr(0,0), tyy.ptr(0,0), &rhsv.ref(0,0));
		_div_dxy_sse<PIN0, PIN1, POUT>(txz.ptr(0,0), tyz.ptr(0,0), &rhsw.ref(0,0));
		_div_dxy_sse<PIN0, PIN1, POUT>(utx.ptr(0,0), uty.ptr(0,0), &rhss.ref(0,0));
	}
	
	void _div_dz(const TempPiZSOAf_ST& tzx0, const TempPiZSOAf_ST& tzy0, const TempPiZSOAf_ST& tzz0, const TempPiZSOAf_ST& utz0,
				 const TempPiZSOAf_ST& tzx1, const TempPiZSOAf_ST& tzy1, const TempPiZSOAf_ST& tzz1, const TempPiZSOAf_ST& utz1)
	{
		static const int PIN0 = TempPiZSOA_ST::PITCH;
		static const int POUT = OutputSOA::PITCH;
		
		_div_dz_sse<PIN0, POUT>(tzx1.ptr(0,0), tzx0.ptr(0,0), &rhsu.ref(0,0));
		_div_dz_sse<PIN0, POUT>(tzy1.ptr(0,0), tzy0.ptr(0,0), &rhsv.ref(0,0));
		_div_dz_sse<PIN0, POUT>(tzz1.ptr(0,0), tzz0.ptr(0,0), &rhsw.ref(0,0));
		_div_dz_sse<PIN0, POUT>(utz1.ptr(0,0), utz0.ptr(0,0), &rhss.ref(0,0));
		}
#endif
	
	void _udot_tx(const InputSOAf_ST& u, const InputSOAf_ST& v, const InputSOAf_ST& w)
	{
#ifndef _SP_COMP_
		printf("DivTensor_QPX::_udot_tx: you should not be here in double precision. Aborting.\n");
		abort();
#else
		static const int PITCHIN = InputSOA_ST::PITCH;
		static const int PITCHTENSOR = TempPiXSOA_ST::PITCH;
		
		const float * const ubase = u.ptr(0,0);
		const float * const vbase = v.ptr(0,0);
		const float * const wbase = w.ptr(0,0);
		
		const float * const txxbase = txx.ptr(0,0);
		const float * const txybase = txy.ptr(0,0);
		const float * const txzbase = txz.ptr(0,0);
		
		float * const utxbase = &utx.ref(0,0);
		
		for(int iy=0; iy<TempPiXSOA_ST::NY; iy++)
		{
			const float * const uptr = ubase + iy*PITCHIN;
			const float * const vptr = vbase + iy*PITCHIN;
			const float * const wptr = wbase + iy*PITCHIN;
			
			const float * const txxptr = txxbase + iy*PITCHTENSOR;
			const float * const txyptr = txybase + iy*PITCHTENSOR;
			const float * const txzptr = txzbase + iy*PITCHTENSOR;
			
			float * const utxptr = utxbase + iy*PITCHTENSOR;
			
			for(int ix=0; ix<TempPiXSOA_ST::NX; ix += 4)
			{
				const vector4double C0 = vec_lda(0L,uptr + ix);
				const vector4double C1 = vec_lda(0L,vptr + ix);
				const vector4double C2 = vec_lda(0L,wptr + ix);

				const vector4double W0 = vec_lda(0L,uptr + ix - 4);
				const vector4double W1 = vec_lda(0L,vptr + ix - 4);
				const vector4double W2 = vec_lda(0L,wptr + ix - 4);

				vec_sta(	 vec_lda(0L,txxptr + ix)*(C0 + LEFT(W0, C0)) +
							 vec_lda(0L,txyptr + ix)*(C1 + LEFT(W1, C1)) +
							 vec_lda(0L,txzptr + ix)*(C2 + LEFT(W2, C2)), 0L, utxptr + ix);
			}
		}
#endif
	}
	
	void _udot_ty(const InputSOAf_ST& u, const InputSOAf_ST& v, const InputSOAf_ST& w)
	{
#ifndef _SP_COMP_
		printf("DivTensor_QPX::_udot_ty: you should not be here in double precision. Aborting.\n");
		abort();
#else
		static const int PITCHIN = InputSOA_ST::PITCH;
		static const int PITCHTENSOR = TempPiYSOA_ST::PITCH;
		
		const float * const ubase = u.ptr(0,0);
		const float * const vbase = v.ptr(0,0);
		const float * const wbase = w.ptr(0,0);
		
		const float * const tyxbase = tyx.ptr(0,0);
		const float * const tyybase = tyy.ptr(0,0);
		const float * const tyzbase = tyz.ptr(0,0);
		
		float * const utybase = &uty.ref(0,0);
		
		for(int iy=0; iy<TempPiYSOA_ST::NY; iy++)
		{
			const float * const uptr = ubase + iy*PITCHIN;
			const float * const vptr = vbase + iy*PITCHIN;
			const float * const wptr = wbase + iy*PITCHIN;
			
			const float * const tyxptr = tyxbase + iy*PITCHTENSOR;
			const float * const tyyptr = tyybase + iy*PITCHTENSOR;
			const float * const tyzptr = tyzbase + iy*PITCHTENSOR;
			
			float * const utyptr = utybase + iy*PITCHTENSOR;
			
			for(int ix=0; ix<TempPiYSOA_ST::NX; ix += 4)
			{
				vec_sta(	 vec_lda(0L,tyxptr + ix)*(vec_lda(0L,uptr + ix) + vec_lda(0L,uptr + ix - PITCHIN)) +
							 vec_lda(0L,tyyptr + ix)*(vec_lda(0L,vptr + ix) + vec_lda(0L,vptr + ix - PITCHIN)) +
							 vec_lda(0L,tyzptr + ix)*(vec_lda(0L,wptr + ix) + vec_lda(0L,wptr + ix - PITCHIN)), 0L, utyptr + ix);
			}
		}
#endif
	}
	
	void _udot_tz(const InputSOAf_ST& u0, const InputSOAf_ST& v0, const InputSOAf_ST& w0,
								 const InputSOAf_ST& u1, const InputSOAf_ST& v1, const InputSOAf_ST& w1,
								 const TempPiZSOAf_ST& tzx, const TempPiZSOAf_ST& tzy, const TempPiZSOAf_ST& tzz,
								 TempPiZSOAf_ST& utz)
	{
		static const int PITCHIN = InputSOA_ST::PITCH;
		static const int PITCHTENSOR = TempPiZSOA_ST::PITCH;
		
		const float * const ubase0 = u0.ptr(0,0);
		const float * const vbase0 = v0.ptr(0,0);
		const float * const wbase0 = w0.ptr(0,0);
		const float * const ubase1 = u1.ptr(0,0);
		const float * const vbase1 = v1.ptr(0,0);
		const float * const wbase1 = w1.ptr(0,0);
		
		const float * const tzxbase = tzx.ptr(0,0);
		const float * const tzybase = tzy.ptr(0,0);
		const float * const tzzbase = tzz.ptr(0,0);
		
		float * const utzbase = &utz.ref(0,0);
		
		for(int iy=0; iy<TempPiZSOA_ST::NY; iy++)
		{
			const float * const uptr0 = ubase0 + iy*PITCHIN;
			const float * const vptr0 = vbase0 + iy*PITCHIN;
			const float * const wptr0 = wbase0 + iy*PITCHIN;
			const float * const uptr1 = ubase1 + iy*PITCHIN;
			const float * const vptr1 = vbase1 + iy*PITCHIN;
			const float * const wptr1 = wbase1 + iy*PITCHIN;
			
			const float * const tzxptr = tzxbase + iy*PITCHTENSOR;
			const float * const tzyptr = tzybase + iy*PITCHTENSOR;
			const float * const tzzptr = tzzbase + iy*PITCHTENSOR;
			
			float * const utzptr = utzbase + iy*PITCHTENSOR;
			
			for(int ix=0; ix<TempPiZSOA_ST::NX; ix += 4)
			{
				vec_sta(	 vec_lda(0L,tzxptr + ix)*(vec_lda(0L,uptr0 + ix) + vec_lda(0L,uptr1 + ix)) +
							 vec_lda(0L,tzyptr + ix)*(vec_lda(0L,vptr0 + ix) + vec_lda(0L,vptr1 + ix)) +
							 vec_lda(0L,tzzptr + ix)*(vec_lda(0L,wptr0 + ix) + vec_lda(0L,wptr1 + ix)), 0L, utzptr + ix);
			}
		}
	}
		
protected:
	
	template <bool accum> inline void write128(float * const dest, const vector4double src)
	{
		if (accum)
			vec_sta(vec_lda(0L,dest) + src, 0L, dest);
		else
			vec_sta(src, 0L, dest);
	}
	
	template<bool accum>
	void _average_xface(const float factor, const TempSOA_ST& src0, const TempSOA_ST& src1, TempPiXSOA_ST& dest)
	{
		static const int PITCHIN = TempSOA_ST::PITCH;
		static const int PITCHOUT = TempPiXSOA_ST::PITCH;
		
		const float * const src0base = src0.ptr(0,0); 
		const float * const src1base = src1.ptr(0,0); 
		float * const destbase = &dest.ref(0,0); 
		
		const vector4double F = vec_splats(factor);
		
		for(int iy=0; iy<TempPiXSOA_ST::NY; iy++)
		{
			const float * const src0ptr = src0base + iy*PITCHIN; 
			const float * const src1ptr = src1base + iy*PITCHIN; 
			float * const destptr = destbase + iy*PITCHOUT;
			
			for(int ix=0; ix<TempPiXSOA_ST::NX; ix+=4)
				write128<accum>(destptr + ix, F*(vec_lda(0L,src0ptr + ix) + vec_lda(0L,src0ptr + ix + PITCHIN) + 
												vec_lda(0L,src1ptr + ix) + vec_lda(0L,src1ptr + ix + PITCHIN)));
		}
	}
	
	template<bool accum>
	void _average_yface(const float factor, const TempSOA_ST& src0, const TempSOA_ST& src1, TempPiYSOA_ST& dest)
	{
		static const int PITCHIN = TempSOA_ST::PITCH;
		static const int PITCHOUT = TempPiYSOA_ST::PITCH;
		
		const float * const src0base = src0.ptr(0,0); 
		const float * const src1base = src1.ptr(0,0); 
		float * const destbase = &dest.ref(0,0); 
		
		const vector4double F = vec_splats(factor);
		
		for(int iy=0; iy<TempPiYSOA_ST::NY; iy++)
		{
			const float * const src0ptr = src0base + iy*PITCHIN; 
			const float * const src1ptr = src1base + iy*PITCHIN; 
			float * const destptr = destbase + iy*PITCHOUT;
						
			vector4double C0 = vec_lda(0L,src0ptr);
			vector4double C1 = vec_lda(0L,src1ptr);
			
			for(int ix=0; ix<TempPiYSOA_ST::NX; ix+=4)
			{
				const vector4double E0 = vec_lda(0L,src0ptr + ix + 4);
				const vector4double E1 = vec_lda(0L,src1ptr + ix + 4);

				write128<accum>(destptr + ix, F*(C0 + RIGHT(C0, E0) + 
												 C1 + RIGHT(C1, E1)));
				C0 = E0;
				C1 = E1;
			}
		}
	}
	
	template<bool accum>
	void _average_zface(const float factor, const TempSOA_ST& src, TempPiZSOA_ST& dest)
	{
		static const int PITCHIN = TempSOA_ST::PITCH;
		static const int PITCHOUT = TempPiZSOA_ST::PITCH;
		
		const float * const srcbase = src.ptr(0,0); 
		float * const destbase = &dest.ref(0,0); 
		
		const vector4double F = vec_splats(factor);
		
		for(int iy=0; iy<TempPiZSOA_ST::NY; iy++)
		{
			const float * const srcptr = srcbase + iy*PITCHIN; 
			float * const destptr = destbase + iy*PITCHOUT;
			
			vector4double C0 = vec_lda(0L,srcptr);
			vector4double C1 = vec_lda(0L,srcptr + PITCHIN);
			
			for(int ix=0; ix<TempPiZSOA_ST::NX; ix+=4)
			{
				const vector4double E0 = vec_lda(0L,srcptr + ix + 4);
				const vector4double E1 = vec_lda(0L,srcptr + ix + PITCHIN + 4);
				
				write128<accum>(destptr + ix, F*(C0 + RIGHT(C0, E0) + 
												 C1 + RIGHT(C1, E1)));
				C0 = E0;
				C1 = E1;
			}
		}
	}
		
public:
	
	DivTensor_QPX(const Real a = 1, const Real dtinvh = 1, const Real h = 1, const Real sigma=1):
	DivTensor_CPP(a, dtinvh, h, sigma)
	{
	}
	
#undef LEFT
#undef RIGHT	

};
