#ifndef _WAVELETSONINTERVALQPX_H_
#define _WAVELETSONINTERVALQPX_H_ 1

#pragma once

#ifdef _SP_COMP_
typedef float Real;
#else
typedef double Real;
#endif

#ifdef _QPXEMU_
#include <xmmintrin.h>
#define vector4double __m128
#define _DIEGO_TRANSPOSE4(a,b,c,d) _MM_TRANSPOSE4_PS(a,b,c,d)
#define vec_lda(a, b) _mm_load_ps(b + ((a)/4))
#define vec_sta(a, b, c) _mm_store_ps(c + ((b)/4), a)
#define vec_mul(a, b) (a) * (b)
#define vec_sub(a, b) (a) - (b)
#define vec_add(a, b) (a) + (b)
#define vec_madd(a, b, c) (a) * (b) + (c)
#define vec_nmsub(a, b, c) (c) - (a) * (b)
#define vec_splats(a) _mm_set1_ps(a)
#endif

#ifndef _DIEGO_TRANSPOSE4
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

namespace WaveletsOnInterval
{
	struct WI4QPX
	{
		template<int NH>
		struct __attribute__((__aligned__(_ALIGNBYTES_))) TinyScratchPad
		{
			FwtAp c[NH][4];
			FwtAp d[NH][4];
			
			void _cpbk(FwtAp *src, FwtAp * dst0, FwtAp * dst1, FwtAp * dst2, FwtAp * dst3)
			{
				for(int is = 0, id = 0; is < NH * 4; is += 16, id += 4)
				{
					vector4double d0 = vec_lda(0, src + is);
					vector4double d1 = vec_lda(sizeof(FwtAp) * 4, src + is);
					vector4double d2 = vec_lda(sizeof(FwtAp) * 4 * 2, src + is);
					vector4double d3 = vec_lda(sizeof(FwtAp) * 4 * 3, src + is);
					
					_DIEGO_TRANSPOSE4(d0, d1, d2, d3);
					
					vec_sta(d0, 0, dst0 + id);
					vec_sta(d1, 0, dst1 + id);
					vec_sta(d2, 0, dst2 + id);
					vec_sta(d3, 0, dst3 + id);
				}
			}
			
			void _put(Real *dst, Real * src0, Real * src1, Real * src2, Real * src3)
			{
				for(int is = 0, id = 0; is < NH; is += 4, id += 16)
				{
					vector4double d0 = vec_lda(0, src0 + is);
					vector4double d1 = vec_lda(0, src1 + is);
					vector4double d2 = vec_lda(0, src2 + is);
					vector4double d3 = vec_lda(0, src3 + is);
					
					_DIEGO_TRANSPOSE4(d0, d1, d2, d3);
					
					vec_sta(d0, 0, dst + id);
					vec_sta(d1, sizeof(Real) * 4, dst + id);
					vec_sta(d2, sizeof(Real) * 4 * 2, dst + id);
					vec_sta(d3, sizeof(Real) * 4 * 3, dst + id);
				}
			}
			
			void copyback(FwtAp * stream0, FwtAp * stream1, FwtAp * stream2, FwtAp * stream3)
			{
				assert(NH % 4 == 0);
				
				_cpbk(&c[0][0], stream0, stream1, stream2, stream3);
				_cpbk(&d[0][0], stream0 + NH, stream1 + NH, stream2 + NH, stream3 + NH);
			}
			
			void put(Real * stream0, Real * stream1, Real * stream2, Real * stream3)
			{
				assert(NH % 4 == 0);
				_put(&c[0][0], stream0, stream1, stream2, stream3);
				_put(&d[0][0], stream0 + NH, stream1 + NH, stream2 + NH, stream3 + NH);
			}
		};
		
		vector4double P_1_16, P_5_16, P_9_16, P_15_16, P_21_16, P_35_16;
		
		WI4QPX()
		{
			P_1_16 = vec_splats(1.f/16);
			P_5_16 = vec_splats(5.f/16);
			P_9_16 = vec_splats(9.f/16);
			P_15_16 = vec_splats(15.f/16);
			P_21_16 = vec_splats(21.f/16);
			P_35_16 = vec_splats(35.f/16);
		}
		
		inline vector4double interp_first(const vector4double f0, const vector4double f1, const vector4double f2, const vector4double f3)
		{
			return vec_madd(P_5_16, f0, vec_madd(P_15_16, f1, vec_nmsub(P_5_16, f2, vec_mul(P_1_16, f3))));
		};
		
		inline vector4double interp_middle(const vector4double f0, const vector4double f1, const vector4double f2, const vector4double f3)
		{
			return vec_nmsub(P_1_16, f0, vec_madd(P_9_16, f1, vec_nmsub(P_1_16, f3, vec_mul(P_9_16, f2))));
		};
		
		inline vector4double interp_onetolast(const vector4double f0, const vector4double f1, const vector4double f2, const vector4double f3)
		{
			return vec_madd(P_1_16, f0, vec_nmsub(P_5_16, f1, vec_madd(P_15_16, f2, vec_mul(P_5_16, f3))));
		};
		
		inline vector4double interp_last(const vector4double f0, const vector4double f1, const vector4double f2, const vector4double f3)
		{
			return vec_nmsub(P_5_16, f0, vec_madd(P_21_16, f1, vec_nmsub(P_35_16, f2, vec_mul(P_35_16, f3))));
		};
		
		template<int N>
		void forward(FwtAp * stream0, FwtAp * stream1, FwtAp * stream2, FwtAp * stream3)
		{		
			enum { NH = N/2 };
			
			assert(N >= 8);
			assert(N % 8 == 0);
			
			TinyScratchPad<NH> mycoeffs;
			
			vector4double d0 = vec_lda(0, stream0);
			vector4double d1 = vec_lda(0, stream1);
			vector4double d2 = vec_lda(0, stream2);
			vector4double d3 = vec_lda(0, stream3);
			
			_DIEGO_TRANSPOSE4(d0, d1, d2, d3);
			
			vector4double d4 = vec_lda(sizeof(FwtAp) * 4, stream0);
			vector4double d5 = vec_lda(sizeof(FwtAp) * 4, stream1);
			vector4double d6 = vec_lda(sizeof(FwtAp) * 4, stream2);
			vector4double d7 = vec_lda(sizeof(FwtAp) * 4, stream3);
			
			_DIEGO_TRANSPOSE4(d4, d5, d6, d7);
			
			vec_sta(d0, 0, mycoeffs.c[0]);
			vec_sta(vec_sub(d1, interp_first(d0, d2, d4, d6)), 0, mycoeffs.d[0]);
			vec_sta(d2, 0, mycoeffs.c[1]);
			vec_sta(vec_sub(d3, interp_middle(d0, d2, d4, d6)), 0, mycoeffs.d[1]);
			
			for(int src = 8, dst = 2; src < N; src += 4, dst += 2)
			{
				vector4double dm2 = d2;
				
				d0 = d4;
				d1 = d5;
				d2 = d6;
				d3 = d7;
				
				d4 = vec_lda(0, stream0 + src);
				d5 = vec_lda(0, stream1 + src);
				d6 = vec_lda(0, stream2 + src);
				d7 = vec_lda(0, stream3 + src);
				
				_DIEGO_TRANSPOSE4(d4, d5, d6, d7);
				
				vec_sta(d0, 0, mycoeffs.c[dst]);
				vec_sta(vec_sub(d1, interp_middle(dm2, d0, d2, d4)), 0, mycoeffs.d[dst]);
				vec_sta(d2, 0, mycoeffs.c[dst + 1]);
				vec_sta(vec_sub(d3, interp_middle(d0, d2, d4, d6)), 0, mycoeffs.d[dst + 1]);
			}
			
			vec_sta(d4, 0, mycoeffs.c[NH-2]);
			vec_sta(vec_sub(d5, interp_onetolast(d0, d2, d4, d6)), 0, mycoeffs.d[NH-2]);
			vec_sta(d6, 0, mycoeffs.c[NH-1]);
			vec_sta(vec_sub(d7, interp_last(d0, d2, d4, d6)), 0, mycoeffs.d[NH-1]);
			
			mycoeffs.copyback(stream0, stream1, stream2, stream3);
		}
		
		template<int N>
		void inverse(Real * stream0, Real * stream1, Real * stream2, Real * stream3)
		{
			enum { NH = N / 2 };
			
			assert(N >= 8);
			assert(N % 8 == 0);
			
			TinyScratchPad<NH> mycoeffs;
			
			mycoeffs.put(stream0, stream1, stream2, stream3);
			
			vector4double c0 = vec_lda(0, mycoeffs.c[0]);
			vector4double c1 = vec_lda(0, mycoeffs.c[1]);
			vector4double c2 = vec_lda(0, mycoeffs.c[2]);
			vector4double c3 = vec_lda(0, mycoeffs.c[3]);
			
			//first part
			{
				vector4double d0 = c0;
				vector4double d1 = vec_add(vec_lda(0, mycoeffs.d[0]), interp_first(c0, c1, c2, c3));
				vector4double d2 = c1;
				vector4double d3 = vec_add(vec_lda(0, mycoeffs.d[1]), interp_middle(c0, c1, c2, c3));
				
				_DIEGO_TRANSPOSE4(d0, d1, d2, d3);
				
				vec_sta(d0, 0, stream0);
				vec_sta(d1, 0, stream1);
				vec_sta(d2, 0, stream2);
				vec_sta(d3, 0, stream3);
			}
			
			for(int dst = 4, src = 2; dst < N - 4; dst += 4, src += 2)
			{
				vector4double cm1 = c1;
				
				c0 = c2;
				c1 = c3;
				
				c2 = vec_lda(0, mycoeffs.c[src + 2]);
				c3 = vec_lda(0, mycoeffs.c[src + 3]);
				
				vector4double d0 = c0;
				vector4double d1 = vec_add(vec_lda(0, mycoeffs.d[src]), interp_middle(cm1, c0, c1, c2));
				vector4double d2 = c1;
				vector4double d3 = vec_add(vec_lda(0, mycoeffs.d[src + 1]), interp_middle(c0, c1, c2, c3));
				
				_DIEGO_TRANSPOSE4(d0, d1, d2, d3);
				
				vec_sta(d0, 0, stream0 + dst);
				vec_sta(d1, 0, stream1 + dst);
				vec_sta(d2, 0, stream2 + dst);
				vec_sta(d3, 0, stream3 + dst);
			}
			
			//last part
			{
				vector4double d0 = c2;
				vector4double d1 = vec_add(vec_lda(0, mycoeffs.d[NH - 2]), interp_onetolast(c0, c1, c2, c3));
				vector4double d2 = c3;
				vector4double d3 = vec_add(vec_lda(0, mycoeffs.d[NH - 1]), interp_last(c0, c1, c2, c3));
				
				_DIEGO_TRANSPOSE4(d0, d1, d2, d3);
				
				vec_sta(d0, 0, stream0 + N - 4);
				vec_sta(d1, 0, stream1 + N - 4);
				vec_sta(d2, 0, stream2 + N - 4);
				vec_sta(d3, 0, stream3 + N - 4);
			}
		}
	};
	
	template<int ROWSIZE, int COLSIZE>
	struct WaveletSweepQPX
	{
		WI4QPX qpxfwt;
		
		template<int BS, bool forward>
		inline void sweep1D(FwtAp data[BS][ROWSIZE])
		{
			if (forward)
				for(int iy = 0; iy < BS; iy += 4)
					qpxfwt.template forward<BS>(&data[iy][0], &data[iy + 1][0], &data[iy + 2][0], &data[iy + 3][0]);
			else 
				for(int iy = 0; iy < BS; iy += 4)
					qpxfwt.template inverse<BS>(&data[iy][0], &data[iy + 1][0], &data[iy + 2][0], &data[iy + 3][0]);
		}
		
		template<int BS>
		inline void _xy_tr_block(int ix, int iy, FwtAp data[BS][ROWSIZE])
		{
			vector4double d0 = vec_lda(0, &data[iy][ix]);
			vector4double d1 = vec_lda(0, &data[iy + 1][ix]);
			vector4double d2 = vec_lda(0, &data[iy + 2][ix]);
			vector4double d3 = vec_lda(0, &data[iy + 3][ix]);
			
			_DIEGO_TRANSPOSE4(d0, d1, d2, d3);
			
			vec_sta(d0, 0, &data[iy][ix]);
			vec_sta(d1, 0, &data[iy + 1][ix]);
			vec_sta(d2, 0, &data[iy + 2][ix]);
			vec_sta(d3, 0, &data[iy + 3][ix]);
		}
		
		template<int BS>
		inline void _xy_tr_2blocks(int ix1, int iy1, int ix2, int iy2, FwtAp data[BS][ROWSIZE])
		{
			vector4double d0 = vec_lda(0, &data[iy1][ix1]);
			vector4double d1 = vec_lda(0, &data[iy1 + 1][ix1]);
			vector4double d2 = vec_lda(0, &data[iy1 + 2][ix1]);
			vector4double d3 = vec_lda(0, &data[iy1 + 3][ix1]);
			
			_DIEGO_TRANSPOSE4(d0, d1, d2, d3);
			
			vector4double d4 = vec_lda(0, &data[iy2][ix2]);
			vector4double d5 = vec_lda(0, &data[iy2 + 1][ix2]);
			vector4double d6 = vec_lda(0, &data[iy2 + 2][ix2]);
			vector4double d7 = vec_lda(0, &data[iy2 + 3][ix2]);
			
			vec_sta(d0, 0, &data[iy2][ix2]);
			vec_sta(d1, 0, &data[iy2 + 1][ix2]);
			vec_sta(d2, 0, &data[iy2 + 2][ix2]);
			vec_sta(d3, 0, &data[iy2 + 3][ix2]);
			
			_DIEGO_TRANSPOSE4(d4, d5, d6, d7);
			
			vec_sta(d4, 0, &data[iy1][ix1]);
			vec_sta(d5, 0, &data[iy1 + 1][ix1]);
			vec_sta(d6, 0, &data[iy1 + 2][ix1]);
			vec_sta(d7, 0, &data[iy1 + 3][ix1]);			
		}
		
		template<int BS>
		inline void xy_transpose(FwtAp data[BS][ROWSIZE])
		{
			for(int iy = 0; iy < BS; iy += 4)
			{
				_xy_tr_block<BS>(iy, iy, data);
				
				for(int ix = iy + 4; ix < BS; ix += 4)
					_xy_tr_2blocks<BS>(ix, iy, iy, ix, data);
			}
		}
		
		template<int BS>
		inline void _xz_tr_block(int ix, int iy, int iz, FwtAp data[BS][COLSIZE][ROWSIZE])
		{
			vector4double d0 = vec_lda(0, &data[iz][iy][ix]);
			vector4double d1 = vec_lda(0, &data[iz + 1][iy][ix]);
			vector4double d2 = vec_lda(0, &data[iz + 2][iy][ix]);
			vector4double d3 = vec_lda(0, &data[iz + 3][iy][ix]);
			
			_DIEGO_TRANSPOSE4(d0, d1, d2, d3);
			
			vec_sta(d0, 0, &data[iz][iy][ix]);
			vec_sta(d1, 0, &data[iz + 1][iy][ix]);
			vec_sta(d2, 0, &data[iz + 2][iy][ix]);
			vec_sta(d3, 0, &data[iz + 3][iy][ix]);
		}
		
		template<int BS>
		inline void _xz_tr_2blocks(int ix1, int iy1, int iz1, int ix2, int iy2, int iz2, FwtAp data[BS][COLSIZE][ROWSIZE])
		{
			vector4double d0 = vec_lda(0, &data[iz1 + 0][iy1][ix1]);
			vector4double d1 = vec_lda(0, &data[iz1 + 1][iy1][ix1]);
			vector4double d2 = vec_lda(0, &data[iz1 + 2][iy1][ix1]);
			vector4double d3 = vec_lda(0, &data[iz1 + 3][iy1][ix1]);
			
			_DIEGO_TRANSPOSE4(d0, d1, d2, d3);
			
			vector4double d4 = vec_lda(0, &data[iz2 + 0][iy2][ix2]);
			vector4double d5 = vec_lda(0, &data[iz2 + 1][iy2][ix2]);
			vector4double d6 = vec_lda(0, &data[iz2 + 2][iy2][ix2]);
			vector4double d7 = vec_lda(0, &data[iz2 + 3][iy2][ix2]);
			
			vec_sta(d0, 0, &data[iz2][iy2][ix2]);
			vec_sta(d1, 0, &data[iz2 + 1][iy2][ix2]);
			vec_sta(d2, 0, &data[iz2 + 2][iy2][ix2]);
			vec_sta(d3, 0, &data[iz2 + 3][iy2][ix2]);
			
			_DIEGO_TRANSPOSE4(d4, d5, d6, d7);
			
			vec_sta(d4, 0, &data[iz1 + 0][iy1][ix1]);
			vec_sta(d5, 0, &data[iz1 + 1][iy1][ix1]);
			vec_sta(d6, 0, &data[iz1 + 2][iy1][ix1]);
			vec_sta(d7, 0, &data[iz1 + 3][iy1][ix1]);			
		}
		
		template<int BS>
		inline void xz_transpose(FwtAp data[BS][COLSIZE][ROWSIZE])
		{
			for(int iy = 0; iy < BS; ++iy)
				for(int iz = 0; iz < BS; iz += 4)
				{
					_xz_tr_block<BS>(iz, iy, iz, data);
					
					for(int ix = iz + 4; ix < BS; ix += 4)
						_xz_tr_2blocks<BS>(ix, iy, iz, iz, iy, ix, data);
				}
		}
		
		template<int BS, bool forward>
		inline void sweep2D(FwtAp data[BS][ROWSIZE])
		{			
			sweep1D<BS, forward>(data);
			xy_transpose<BS>(data);
			sweep1D<BS, forward>(data);
		}
		
		template<int BS, bool bForward>
		inline void sweep3D(FwtAp data[BS][COLSIZE][ROWSIZE])
		{			
			if(bForward)
			{
				for(int iz = 0; iz < BS; ++iz)
					sweep2D<BS, true>(data[iz]);
				
				xz_transpose<BS>(data);
				
				for(int iz = 0; iz < BS; ++iz)
					for(int iy = 0; iy < BS; iy += 4)
						qpxfwt.template forward<BS>(&data[iz][iy][0], &data[iz][iy + 1][0], &data[iz][iy + 2][0], &data[iz][iy + 3][0]);
			}
			else
			{			
				for(int iz = 0; iz < BS; ++iz)
					for(int iy = 0; iy < BS; iy += 4)
						qpxfwt.template inverse<BS>(&data[iz][iy][0], &data[iz][iy + 1][0], &data[iz][iy + 2][0], &data[iz][iy + 3][0]);
				
				xz_transpose<BS>(data);
				
				for(int iz = 0; iz < BS; ++iz)
					sweep2D<BS, false>(data[iz]);
			}
		}
	};
}

/* I USE THIS TO TEST MY QPX EMULATION BECAUSE IBM DOES NOT GIVE ME ACCESS TO BGQ :(
#include "WaveletsOnInterval.h"
#include "Timer.h"
struct TestFWTQPX
{
	TestFWTQPX()
	{
		enum {N = 32 } ;
		printf("heeeelllllllloooooooooo\n"); 
		WaveletsOnInterval::WI4QPX asdo;
		WaveletsOnInterval::WaveletSweepQPX<N, N> mysweep;
		WaveletsOnInterval::FwtAp __attribute__((aligned(16))) input[4][N];
		WaveletsOnInterval::FwtAp  __attribute__((aligned(16))) output[4][N];
		
		WaveletsOnInterval::FwtAp __attribute__((aligned(16))) inputM[N][N];
		WaveletsOnInterval::FwtAp __attribute__((aligned(16))) outputM[N][N];

		WaveletsOnInterval::FwtAp __attribute__((aligned(16))) inputM3D[N][N][N];
		WaveletsOnInterval::FwtAp __attribute__((aligned(16))) outputM3D[N][N][N];
		
		for(int d = 0; d < 4; ++d)
			for(int i = 0; i < N; ++i)
				input[d][i] = output[d][i] = 1 + drand48();
		
		for(int d = 0; d < N; ++d)
			for(int i = 0; i < N; ++i)
				inputM[d][i] = outputM[d][i] = 1 + drand48();

		for(int k = 0; k < N; ++k)		
		for(int d = 0; d < N; ++d)
			for(int i = 0; i < N; ++i)
				inputM3D[k][d][i] = outputM3D[k][d][i] = 1 + drand48();
		
		mysweep.xy_transpose<N>(outputM);
		
		for(int d = 0; d < N; ++d)
			for(int i = 0; i < N; ++i)
			{
			//	printf("%d %d\n", d, i);
				assert(inputM[i][d] == outputM[d][i]);
			}
		
		mysweep.xz_transpose<N>(outputM3D);

		for(int k = 0; k < N; ++k)		
			for(int d = 0; d < N; ++d)
				for(int i = 0; i < N; ++i)
					assert(inputM3D[i][k][d] == outputM3D[d][k][i]);

//		exit(0);
	
		WaveletsOnInterval::WI4<false> myref;
		Timer timer; timer.start();
	
		const int NTIMES = 1e4;
		for(int i = 0; i< NTIMES; ++i)
		{
#if 0
			for(int d = 0; d < 4; ++d)
			{
				myref.transform<N, true>(output[d]);
				myref.transform<N, false>(output[d]);
			}
#elif 0
			asdo.forward<N>(output[0], output[1], output[2], output[3]);
			asdo.inverse<N>(output[0], output[1], output[2], output[3]);
#else
			//mysweep.xy_transpose<N>(outputM);
			mysweep.xz_transpose<N>(outputM3D);
#endif
		}
				
		printf("spent: %f s\n", timer.stop());
		
		for(int d = 0; d < 4; ++d)
			for(int i = 0; i < N; ++i)
			{
				if (fabs(input[d][i]- output[d][i])>=1e-6) printf("error: %.2e -> input[%d %d] = %f output = %f\n", input[d][i]- output[d][i], d, i, input[d][i], output[d][i]);
			}
		
		exit(0);
	}
};

static TestFWTQPX culo;
 */

#endif
