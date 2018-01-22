/*
 *  DivSOA2D_QPX.cpp
 *  MPCFcore
 *
 *  Created by Fabian Wermelinger 03/03/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#include "DivSOA2D_QPX.h"

#define MYCODE vec_gpci(01234)
#define SEMIPOSITIVE_K

inline void _qpx_xrhsadd(Real * const f, Real * const r)
{
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix+=4)
		{
			Real * const entry_in = f + ix + TempSOA::PITCH*iy;

			const vector4double data0 = vec_lda(0L, entry_in);
			const vector4double data1 = vec_lda(sizeof(Real) * 4, entry_in);
			const vector4double datanext = vec_perm(data0, data1, MYCODE);

			vec_sta(vec_sub(datanext, data0), 0L, r + ix + OutputSOA::PITCH*iy);
		}
}

inline void _qpx_yrhsadd(Real * const f, Real * const r)
{
	enum {
		SP = TempSOA::PITCH,
		DP = OutputSOA::PITCH
	} ;

	for(int iy=0; iy<OutputSOA::NY; iy+=4)
		for(int ix=0; ix<OutputSOA::NX; ix+=4)
		{
			const int offset_in = iy + SP * ix;
			const int offset_out = ix + DP * iy;

			const vector4double dataA0 = vec_lda(0L, f + offset_in);
			const vector4double dataA1 = vec_lda(sizeof(Real) * 4, f + offset_in);
			const vector4double dataB0 = vec_lda(0L, f + offset_in + SP);
			const vector4double dataB1 = vec_lda(sizeof(Real) * 4, f + offset_in + SP);
			const vector4double dataC0 = vec_lda(0L, f + offset_in + 2 * SP);
			const vector4double dataC1 = vec_lda(sizeof(Real) * 4, f + offset_in + 2 * SP);
			const vector4double dataD0 = vec_lda(0L, f + offset_in + 3 * SP);
			const vector4double dataD1 = vec_lda(sizeof(Real) * 4, f + offset_in + 3 * SP);

			vector4double rhs0, rhs1, rhs2, rhs3;

			rhs0 = vec_sub(vec_perm(dataA0, dataA1, MYCODE), dataA0);
			rhs1 = vec_sub(vec_perm(dataB0, dataB1, MYCODE), dataB0);
			rhs2 = vec_sub(vec_perm(dataC0, dataC1, MYCODE), dataC0);
			rhs3 = vec_sub(vec_perm(dataD0, dataD1, MYCODE), dataD0);

			_DIEGO_TRANSPOSE4(rhs0, rhs1, rhs2, rhs3);

			rhs0 = vec_add(rhs0, vec_lda(0L, r + offset_out));
			rhs1 = vec_add(rhs1, vec_lda(0L, r + offset_out + DP));
			rhs2 = vec_add(rhs2, vec_lda(0L, r + offset_out + 2 * DP));
			rhs3 = vec_add(rhs3, vec_lda(0L, r + offset_out + 3 * DP));

			vec_sta(rhs0, 0L, r + offset_out);
			vec_sta(rhs1, 0L, r + offset_out + DP);
			vec_sta(rhs2, 0L, r + offset_out + 2 * DP);
			vec_sta(rhs3, 0L, r + offset_out + 3 * DP);
		}
}
inline void _qpx_zrhsadd(Real * const fb, Real * const ff, Real * const r)
{
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix+=4)
		{
			const vector4double result = vec_sub(vec_lda(0L, ff + ix + TempSOA::PITCH*iy),
                                                 vec_lda(0L, fb + ix + TempSOA::PITCH*iy));

			const vector4double oldvalue = vec_lda(0L, r + ix + OutputSOA::PITCH*iy);

			vec_sta(vec_add(result, oldvalue), 0L, r + ix + OutputSOA::PITCH*iy);
		}
}

inline void _qpx_xextraterm(Real * const xm, Real * const xp, Real * const output)
{
	for(int iy=0; iy<OutputSOA::NY; ++iy)
	{
		Real * const xmptr = xm + iy * TempSOA::PITCH;
		Real * const xpptr = xp + iy * TempSOA::PITCH;
		Real * const outputptr = output + iy * OutputSOA::PITCH;

		for(int ix=0; ix<OutputSOA::NX; ix += 4)
		{
			const vector4double datam = vec_perm(vec_lda(0L, xmptr + ix),
                                                 vec_lda(sizeof(Real) * 4, xmptr + ix),
                                                 vec_gpci(01234));

			const vector4double datap = vec_lda(0L, xpptr + ix);

			vec_sta(vec_add(datam, datap), 0L, outputptr + ix);
		}
	}
}


inline void _qpx_xextraterm(Real * const _um, Real * const _up, Real * const _am, Real * const _ap, Real * const _as, Real * const output)
{
    const vector4double F_0 = vec_splats((Real)0.0);
    const vector4double F_1 = vec_splats((Real)1.0);
    const vector4double F_m1 = vec_splats((Real)-1.0);
    const vector4double F_1_2 = vec_splats((Real)0.5);
    const vector4double F_m1_2 = vec_splats((Real)-0.5);

	for(int iy=0; iy<OutputSOA::NY; ++iy)
	{
		Real * const amptr = _am + iy * TempSOA::PITCH;
		Real * const apptr = _ap + iy * TempSOA::PITCH;
		Real * const asptr = _as + iy * TempSOA::PITCH;
        Real * const umptr = _um + iy * TempSOA::PITCH;
		Real * const upptr = _up + iy * TempSOA::PITCH;

		Real * const outputptr = output + iy * OutputSOA::PITCH;

		for(int ix=0; ix<OutputSOA::NX; ix += 4)
		{
			const vector4double am0 = vec_lda(0L, amptr + ix);
			const vector4double am1 = vec_perm(am0, vec_lda(4 * sizeof(Real), amptr + ix), MYCODE);
			const vector4double ap0 = vec_lda(0L, apptr + ix);
			const vector4double ap1 = vec_perm(ap0, vec_lda(4 * sizeof(Real), apptr + ix), MYCODE);
            const vector4double as0 = vec_lda(0L, asptr + ix);
            const vector4double as1 = vec_perm(as0, vec_lda(4 * sizeof(Real), asptr + ix), MYCODE);

			const vector4double um0 = vec_lda(0L, umptr + ix);
			const vector4double um1 = vec_perm(um0, vec_lda(4 * sizeof(Real), umptr + ix), MYCODE);
			const vector4double up0 = vec_lda(0L, upptr + ix);
			const vector4double up1 = vec_perm(up0, vec_lda(4 * sizeof(Real), upptr + ix), MYCODE);

            vector4double sign_star = vec_sel(vec_cmpgt(as1,F_0), F_0, vec_cmpeq(vec_abs(as1),F_0));
            /* vector4double sign_star = vec_sel(F_m1,F_1,vec_cmpgt(as1, F_0)); */
            vector4double s_minus = mymin(am1, F_0);
            vector4double s_pluss = mymax(ap1, F_0);
            vector4double chi_starm = vec_sub(vec_mul(vec_sub(am1,um1),myreciprocal<preclevel>(vec_sub(am1,as1))), F_1);
            vector4double chi_starp = vec_sub(vec_mul(vec_sub(ap1,up1),myreciprocal<preclevel>(vec_sub(ap1,as1))), F_1);

            vector4double term2 = vec_mul(vec_madd(F_m1_2,sign_star,F_1_2),vec_madd(s_pluss, chi_starp, up1));
            const vector4double u_hllc1 = vec_madd(vec_madd(F_1_2,sign_star,F_1_2), vec_madd(s_minus,chi_starm,um1), term2);

            sign_star = vec_sel(vec_cmpgt(as0,F_0), F_0, vec_cmpeq(vec_abs(as0),F_0));
            /* sign_star = vec_sel(F_m1,F_1,vec_cmpgt(as0, F_0)); */
            s_minus = mymin(am0, F_0);
            s_pluss = mymax(ap0, F_0);
            chi_starm = vec_sub(vec_mul(vec_sub(am0,um0),myreciprocal<preclevel>(vec_sub(am0,as0))), F_1);
            chi_starp = vec_sub(vec_mul(vec_sub(ap0,up0),myreciprocal<preclevel>(vec_sub(ap0,as0))), F_1);

            term2 = vec_mul(vec_madd(F_m1_2,sign_star,F_1_2),vec_madd(s_pluss, chi_starp, up0));
            const vector4double u_hllc0 = vec_madd(vec_madd(F_1_2,sign_star,F_1_2), vec_madd(s_minus,chi_starm,um0), term2);

			const vector4double result = vec_sub(u_hllc1,u_hllc0);

			vec_sta(result, 0L, outputptr + ix);
		}
	}
}

inline void _qpx_yextraterm(Real * const xm, Real * const xp, Real * const output)
{
	enum {
		SP = TempSOA::PITCH,
		DP = OutputSOA::PITCH,
		JUMP = sizeof(Real) * 4
	};

	for(int iy=0; iy<OutputSOA::NY; iy+=4)
		for(int ix=0; ix<OutputSOA::NX; ix+=4)
		{
			const int offset_in0 = iy + SP * ix;
			const int offset_in1 = iy + SP * (ix + 1);
			const int offset_in2 = iy + SP * (ix + 2);
			const int offset_in3 = iy + SP * (ix + 3);

			const int offset_out = ix + DP * iy;

			vector4double rhs0 = vec_add(vec_lda(0L, xp + offset_in0), vec_perm(vec_lda(0L, xm + offset_in0), vec_lda(JUMP, xm + offset_in0), MYCODE));
			vector4double rhs1 = vec_add(vec_lda(0L, xp + offset_in1), vec_perm(vec_lda(0L, xm + offset_in1), vec_lda(JUMP, xm + offset_in1), MYCODE));
			vector4double rhs2 = vec_add(vec_lda(0L, xp + offset_in2), vec_perm(vec_lda(0L, xm + offset_in2), vec_lda(JUMP, xm + offset_in2), MYCODE));
			vector4double rhs3 = vec_add(vec_lda(0L, xp + offset_in3), vec_perm(vec_lda(0L, xm + offset_in3), vec_lda(JUMP, xm + offset_in3), MYCODE));

			_DIEGO_TRANSPOSE4(rhs0, rhs1, rhs2, rhs3);

			rhs0 = vec_add(rhs0, vec_lda(0L, output + offset_out));
			rhs1 = vec_add(rhs1, vec_lda(0L, output + offset_out + DP));
			rhs2 = vec_add(rhs2, vec_lda(0L, output + offset_out + 2 * DP));
			rhs3 = vec_add(rhs3, vec_lda(0L, output + offset_out + 3 * DP));

			vec_sta(rhs0, 0L, output + offset_out);
			vec_sta(rhs1, 0L, output + offset_out + DP);
			vec_sta(rhs2, 0L, output + offset_out + 2 * DP);
			vec_sta(rhs3, 0L, output + offset_out + 3 * DP);
		}
}

struct __attribute__((__aligned__(_ALIGNBYTES_)))  TinyScratchPadExtraTerm { Real tmp[4][4];};

inline void _qpx_yextraterm(Real * const _um, Real * const _up, Real * const _am, Real * const _ap, Real * const _as, Real * const output)
{
    const vector4double F_0 = vec_splats((Real)0.0);
    const vector4double F_1 = vec_splats((Real)1.0);
    const vector4double F_m1 = vec_splats((Real)-1.0);
    const vector4double F_1_2 = vec_splats((Real)0.5);
    const vector4double F_m1_2 = vec_splats((Real)-0.5);

	enum {
		SP = TempSOA::PITCH,
		DP = OutputSOA::PITCH,
		JUMP = sizeof(Real) * 4
	};

	TinyScratchPadExtraTerm scratchpad;

	Real * const tmp = (Real *)&scratchpad.tmp[0][0];

	for(int iy=0; iy<OutputSOA::NY; iy+=4)
	{
		for(int ix=0; ix<OutputSOA::NX; ix+=4)
		{
			for(int j=0; j<4; ++j)
			{
				const int offset_in = iy + SP * (ix + j);

				const vector4double am0 = vec_lda(0L, _am + offset_in);
				const vector4double am1 = vec_perm(am0, vec_lda(JUMP, _am + offset_in), MYCODE);
				const vector4double ap0 = vec_lda(0L, _ap + offset_in);
				const vector4double ap1 = vec_perm(ap0, vec_lda(JUMP, _ap + offset_in), MYCODE);
                const vector4double as0 = vec_lda(0L, _as + offset_in);
				const vector4double as1 = vec_perm(as0, vec_lda(JUMP, _as + offset_in), MYCODE);

				const vector4double um0 = vec_lda(0L, _um + offset_in);
				const vector4double um1 = vec_perm(um0, vec_lda(JUMP, _um + offset_in), MYCODE);
				const vector4double up0 = vec_lda(0L, _up + offset_in);
				const vector4double up1 = vec_perm(up0, vec_lda(JUMP, _up + offset_in), MYCODE);

                vector4double sign_star = vec_sel(vec_cmpgt(as1,F_0), F_0, vec_cmpeq(vec_abs(as1),F_0));
                /* vector4double sign_star = vec_sel(F_m1,F_1,vec_cmpgt(as1, F_0)); */
                vector4double s_minus = mymin(am1, F_0);
                vector4double s_pluss = mymax(ap1, F_0);
                vector4double chi_starm = vec_sub(vec_mul(vec_sub(am1,um1),myreciprocal<preclevel>(vec_sub(am1,as1))), F_1);
                vector4double chi_starp = vec_sub(vec_mul(vec_sub(ap1,up1),myreciprocal<preclevel>(vec_sub(ap1,as1))), F_1);

                vector4double term2 = vec_mul(vec_madd(F_m1_2,sign_star,F_1_2),vec_madd(s_pluss, chi_starp, up1));
                const vector4double u_hllc1 = vec_madd(vec_madd(F_1_2,sign_star,F_1_2), vec_madd(s_minus,chi_starm,um1), term2);

                sign_star = vec_sel(vec_cmpgt(as0,F_0), F_0, vec_cmpeq(vec_abs(as0),F_0));
                /* sign_star = vec_sel(F_m1,F_1,vec_cmpgt(as0, F_0)); */
                s_minus = mymin(am0, F_0);
                s_pluss = mymax(ap0, F_0);
                chi_starm = vec_sub(vec_mul(vec_sub(am0,um0),myreciprocal<preclevel>(vec_sub(am0,as0))), F_1);
                chi_starp = vec_sub(vec_mul(vec_sub(ap0,up0),myreciprocal<preclevel>(vec_sub(ap0,as0))), F_1);

                term2 = vec_mul(vec_madd(F_m1_2,sign_star,F_1_2),vec_madd(s_pluss, chi_starp, up0));
                const vector4double u_hllc0 = vec_madd(vec_madd(F_1_2,sign_star,F_1_2), vec_madd(s_minus,chi_starm,um0), term2);

                const vector4double result = vec_sub(u_hllc1,u_hllc0);

				vec_sta(result, 0L, tmp + 4 * j);
			}

			{
				vector4double sum0 = vec_lda(0L, tmp);
				vector4double sum1 = vec_lda(JUMP, tmp);
				vector4double sum2 = vec_lda(2 * JUMP, tmp);
				vector4double sum3 = vec_lda(3 * JUMP, tmp);

				_DIEGO_TRANSPOSE4(sum0, sum1, sum2, sum3);

				const int offset_out = ix + DP * iy;

				sum0 = vec_add(sum0, vec_lda(0L, output + offset_out));
				sum1 = vec_add(sum1, vec_lda(0L, output + offset_out + DP));
				sum2 = vec_add(sum2, vec_lda(0L, output + offset_out + 2 * DP));
				sum3 = vec_add(sum3, vec_lda(0L, output + offset_out + 3 * DP));

				vec_sta(sum0, 0L, output + offset_out);
				vec_sta(sum1, 0L, output + offset_out + DP);
				vec_sta(sum2, 0L, output + offset_out + 2 * DP);
				vec_sta(sum3, 0L, output + offset_out + 3 * DP);
			}
		}
	}
}

inline void _qpx_zextraterm(Real * const xm, Real * const xp, Real * const output)
{
	for(int iy=0; iy<OutputSOA::NY; ++iy)
	{
		Real * const xmptr = xm + iy * TempSOA::PITCH;
		Real * const xpptr = xp + iy * TempSOA::PITCH;
		Real * const outputptr = output + iy * OutputSOA::PITCH;

		for(int ix=0; ix<OutputSOA::NX; ix += 4)
		{
			const vector4double datam = vec_lda(0L, xmptr + ix);
			const vector4double datap = vec_lda(0L, xpptr + ix);
			const vector4double oldvalue = vec_lda(0L, outputptr + ix);

			vec_sta(vec_add(oldvalue, vec_add(datam, datap)), 0L, outputptr + ix);
		}
	}
}

inline void _qpx_zextraterm(Real * const _um0, Real * const _up0, Real * const _um1, Real * const _up1,
                            Real * const _am0, Real * const _ap0, Real * const _am1, Real * const _ap1,
                            Real * const _as0, Real * const _as1,
                            Real * const output)
{
    const vector4double F_0 = vec_splats((Real)0.0);
    const vector4double F_1 = vec_splats((Real)1.0);
    const vector4double F_m1 = vec_splats((Real)-1.0);
    const vector4double F_1_2 = vec_splats((Real)0.5);
    const vector4double F_m1_2 = vec_splats((Real)-0.5);

	for(int iy=0; iy<OutputSOA::NY; ++iy)
	{
		Real * const amptr0 = _am0 + iy * TempSOA::PITCH;
		Real * const apptr0 = _ap0 + iy * TempSOA::PITCH;
		Real * const asptr0 = _as0 + iy * TempSOA::PITCH;
		Real * const umptr0 = _um0 + iy * TempSOA::PITCH;
		Real * const upptr0 = _up0 + iy * TempSOA::PITCH;

		Real * const amptr1 = _am1 + iy * TempSOA::PITCH;
		Real * const apptr1 = _ap1 + iy * TempSOA::PITCH;
		Real * const asptr1 = _as1 + iy * TempSOA::PITCH;
		Real * const umptr1 = _um1 + iy * TempSOA::PITCH;
		Real * const upptr1 = _up1 + iy * TempSOA::PITCH;

		Real * const outputptr = output + iy * OutputSOA::PITCH;

		for(int ix=0; ix<OutputSOA::NX; ix += 4)
		{
			const vector4double am0 = vec_lda(0L, amptr0 + ix);
			const vector4double ap0 = vec_lda(0L, apptr0 + ix);
			const vector4double as0 = vec_lda(0L, asptr0 + ix);
			const vector4double um0 = vec_lda(0L, umptr0 + ix);
			const vector4double up0 = vec_lda(0L, upptr0 + ix);

			const vector4double am1 = vec_lda(0L, amptr1 + ix);
			const vector4double ap1 = vec_lda(0L, apptr1 + ix);
			const vector4double as1 = vec_lda(0L, asptr1 + ix);
			const vector4double um1 = vec_lda(0L, umptr1 + ix);
			const vector4double up1 = vec_lda(0L, upptr1 + ix);

            vector4double sign_star = vec_sel(vec_cmpgt(as1,F_0), F_0, vec_cmpeq(vec_abs(as1),F_0));
            /* vector4double sign_star = vec_sel(F_m1,F_1,vec_cmpgt(as1, F_0)); */
            vector4double s_minus = mymin(am1, F_0);
            vector4double s_pluss = mymax(ap1, F_0);
            vector4double chi_starm = vec_sub(vec_mul(vec_sub(am1,um1),myreciprocal<preclevel>(vec_sub(am1,as1))), F_1);
            vector4double chi_starp = vec_sub(vec_mul(vec_sub(ap1,up1),myreciprocal<preclevel>(vec_sub(ap1,as1))), F_1);

            vector4double term2 = vec_mul(vec_madd(F_m1_2,sign_star,F_1_2),vec_madd(s_pluss, chi_starp, up1));
            const vector4double u_hllc1 = vec_madd(vec_madd(F_1_2,sign_star,F_1_2), vec_madd(s_minus,chi_starm,um1), term2);

            sign_star = vec_sel(vec_cmpgt(as0,F_0), F_0, vec_cmpeq(vec_abs(as0),F_0));
            /* sign_star = vec_sel(F_m1,F_1,vec_cmpgt(as0, F_0)); */
            s_minus = mymin(am0, F_0);
            s_pluss = mymax(ap0, F_0);
            chi_starm = vec_sub(vec_mul(vec_sub(am0,um0),myreciprocal<preclevel>(vec_sub(am0,as0))), F_1);
            chi_starp = vec_sub(vec_mul(vec_sub(ap0,up0),myreciprocal<preclevel>(vec_sub(ap0,as0))), F_1);

            term2 = vec_mul(vec_madd(F_m1_2,sign_star,F_1_2),vec_madd(s_pluss, chi_starp, up0));
            const vector4double u_hllc0 = vec_madd(vec_madd(F_1_2,sign_star,F_1_2), vec_madd(s_minus,chi_starm,um0), term2);

            const vector4double result = vec_sub(u_hllc1,u_hllc0);
			const vector4double oldvalue = vec_lda(0L, outputptr + ix);

			vec_sta(vec_add(oldvalue, result), 0L, outputptr + ix);
		}
	}
}

inline vector4double _computeK_lambda(const vector4double alpha2, const vector4double pressure,
        const vector4double gamma1, const vector4double gamma2,
        const vector4double pc1, const vector4double pc2)
{
    const vector4double gpc1 = vec_mul(gamma1, vec_add(pressure, pc1));
    const vector4double gpc2 = vec_mul(gamma2, vec_add(pressure, pc2));

    const vector4double sela2   = vec_cmpgt(alpha2, vec_splats(0.0f));
    const vector4double selgpc1 = vec_cmpgt(gpc1,   vec_splats(0.0f));
    const vector4double selgpc2 = vec_cmpgt(gpc2,   vec_splats(0.0f));

    const vector4double a2_pos   = vec_sel(vec_splats(1.0f), alpha2, sela2);
    const vector4double gpc1_pos = vec_sel(vec_splats(1.0f), gpc1,   selgpc1);
    const vector4double gpc2_pos = vec_sel(vec_splats(0.0f), gpc2,   selgpc2);
    const vector4double a1_pos   = vec_sub(vec_splats(1.0f), a2_pos);
    const vector4double lambda   = vec_mul(vec_madd(vec_splats(0.5f),sela2,vec_splats(0.5f)), vec_madd(vec_splats(0.5f),selgpc1,vec_splats(0.5f)));

    return vec_mul(vec_mul(vec_mul(vec_mul(lambda,a1_pos),a2_pos), vec_sub(gpc1_pos, gpc2_pos)), myreciprocal<preclevel>(vec_add(vec_mul(a1_pos,gpc2_pos), vec_mul(a2_pos,gpc1_pos))));
}

inline void _qpx_xextraterm_K(Real * const a2m, Real * const a2p, Real * const pm, Real * const pp,  Real * const output,
        const Real g1, const Real g2, const Real pc1, const Real pc2)
{
    const vector4double F_1 = vec_splats(1.0f);
    const vector4double G1  = vec_splats(g1);
    const vector4double G2  = vec_splats(g2);
    const vector4double P1  = vec_splats(pc1);
    const vector4double P2  = vec_splats(pc2);

    for(int iy=0; iy<OutputSOA::NY; ++iy)
    {
        Real * const a2mptr = a2m + iy * TempSOA::PITCH;
        Real * const a2pptr = a2p + iy * TempSOA::PITCH;
        Real * const pmptr = pm + iy * TempSOA::PITCH;
        Real * const ppptr = pp + iy * TempSOA::PITCH;
        Real * const outputptr = output + iy * OutputSOA::PITCH;

        for(int ix=0; ix<OutputSOA::NX; ix += 4)
        {
            const vector4double a2m0 = vec_perm(vec_lda(0L, a2mptr + ix), vec_lda(sizeof(Real) * 4, a2mptr + ix), MYCODE);
            const vector4double pm0  = vec_perm(vec_lda(0L, pmptr + ix), vec_lda(sizeof(Real) * 4, pmptr + ix), MYCODE);
            const vector4double a2p0 = vec_lda(0L, a2pptr + ix);
            const vector4double pp0  = vec_lda(0L, ppptr + ix);

#ifndef SEMIPOSITIVE_K /* Does NOT check for semi-positivity of \alpha_k and \rho_k\c_k^2 */
            const vector4double a1m0 = vec_sub(F_1, a2m0);
            const vector4double Km0  = vec_mul(vec_mul(vec_mul(a1m0,a2m0), vec_sub(vec_mul(G1,vec_add(pm0,P1)),vec_mul(G2,vec_add(pm0,P2)))), myreciprocal<preclevel>(vec_add(vec_mul(vec_mul(a1m0,G2),vec_add(pm0,P2)), vec_mul(vec_mul(a2m0,G1),vec_add(pm0,P1)))));

            const vector4double a1p0 = vec_sub(F_1, a2p0);
            const vector4double Kp0  = vec_mul(vec_mul(vec_mul(a1p0,a2p0), vec_sub(vec_mul(G1,vec_add(pp0,P1)),vec_mul(G2,vec_add(pp0,P2)))), myreciprocal<preclevel>(vec_add(vec_mul(vec_mul(a1p0,G2),vec_add(pp0,P2)), vec_mul(vec_mul(a2p0,G1),vec_add(pp0,P1)))));
#else /* Does check for semi-positivity of \alpha_k and \rho_k\c_k^2 */
            const vector4double Km0  = _computeK_lambda(a2m0, pm0, G1, G2, P1, P2); // i+1
            const vector4double Kp0  = _computeK_lambda(a2p0, pp0, G1, G2, P1, P2); // i
#endif /* SEMIPOSITIVE_K */
            vec_sta(vec_add(Kp0, Km0), 0L, outputptr + ix);
        }
    }
}

inline void _qpx_yextraterm_K(Real * const a2m, Real * const a2p, Real * const pm, Real * const pp,  Real * const output,
        const Real g1, const Real g2, const Real pc1, const Real pc2)
{
    enum {
        SP = TempSOA::PITCH,
        DP = OutputSOA::PITCH,
        JUMP = sizeof(Real) * 4
    };

    const vector4double F_1 = vec_splats(1.0f);
    const vector4double G1  = vec_splats(g1);
    const vector4double G2  = vec_splats(g2);
    const vector4double P1  = vec_splats(pc1);
    const vector4double P2  = vec_splats(pc2);

    for(int iy=0; iy<OutputSOA::NY; iy+=4)
        for(int ix=0; ix<OutputSOA::NX; ix+=4)
        {
            const int offset_in0 = iy + SP * ix;
            const int offset_in1 = iy + SP * (ix + 1);
            const int offset_in2 = iy + SP * (ix + 2);
            const int offset_in3 = iy + SP * (ix + 3);

            const int offset_out = ix + DP * iy;

            // load registers
            // column 0
            const vector4double a2m0 = vec_perm(vec_lda(0L, a2m + offset_in0), vec_lda(JUMP, a2m + offset_in0), MYCODE);
            const vector4double pm0  = vec_perm(vec_lda(0L, pm + offset_in0), vec_lda(JUMP, pm + offset_in0), MYCODE);
            const vector4double a2p0 = vec_lda(0L, a2p + offset_in0);
            const vector4double pp0  = vec_lda(0L, pp + offset_in0);
            // column 1
            const vector4double a2m1 = vec_perm(vec_lda(0L, a2m + offset_in1), vec_lda(JUMP, a2m + offset_in1), MYCODE);
            const vector4double pm1  = vec_perm(vec_lda(0L, pm + offset_in1), vec_lda(JUMP, pm + offset_in1), MYCODE);
            const vector4double a2p1 = vec_lda(0L, a2p + offset_in1);
            const vector4double pp1  = vec_lda(0L, pp + offset_in1);
            // column 2
            const vector4double a2m2 = vec_perm(vec_lda(0L, a2m + offset_in2), vec_lda(JUMP, a2m + offset_in2), MYCODE);
            const vector4double pm2  = vec_perm(vec_lda(0L, pm + offset_in2), vec_lda(JUMP, pm + offset_in2), MYCODE);
            const vector4double a2p2 = vec_lda(0L, a2p + offset_in2);
            const vector4double pp2  = vec_lda(0L, pp + offset_in2);
            // column 3
            const vector4double a2m3 = vec_perm(vec_lda(0L, a2m + offset_in3), vec_lda(JUMP, a2m + offset_in3), MYCODE);
            const vector4double pm3  = vec_perm(vec_lda(0L, pm + offset_in3), vec_lda(JUMP, pm + offset_in3), MYCODE);
            const vector4double a2p3 = vec_lda(0L, a2p + offset_in3);
            const vector4double pp3  = vec_lda(0L, pp + offset_in3);

#ifndef SEMIPOSITIVE_K /* Does NOT check for semi-positivity of \alpha_k and \rho_k\c_k^2 */
            // column 0
            const vector4double a1m0 = vec_sub(F_1, a2m0);
            const vector4double Km0  = vec_mul(vec_mul(vec_mul(a1m0,a2m0), vec_sub(vec_mul(G1,vec_add(pm0,P1)),vec_mul(G2,vec_add(pm0,P2)))), myreciprocal<preclevel>(vec_add(vec_mul(vec_mul(a1m0,G2),vec_add(pm0,P2)), vec_mul(vec_mul(a2m0,G1),vec_add(pm0,P1)))));

            const vector4double a1p0 = vec_sub(F_1, a2p0);
            const vector4double Kp0  = vec_mul(vec_mul(vec_mul(a1p0,a2p0), vec_sub(vec_mul(G1,vec_add(pp0,P1)),vec_mul(G2,vec_add(pp0,P2)))), myreciprocal<preclevel>(vec_add(vec_mul(vec_mul(a1p0,G2),vec_add(pp0,P2)), vec_mul(vec_mul(a2p0,G1),vec_add(pp0,P1)))));

            vector4double rhs0 = vec_add(Kp0, Km0);

            // column 1
            const vector4double a1m1 = vec_sub(F_1, a2m1);
            const vector4double Km1  = vec_mul(vec_mul(vec_mul(a1m1,a2m1), vec_sub(vec_mul(G1,vec_add(pm1,P1)),vec_mul(G2,vec_add(pm1,P2)))), myreciprocal<preclevel>(vec_add(vec_mul(vec_mul(a1m1,G2),vec_add(pm1,P2)), vec_mul(vec_mul(a2m1,G1),vec_add(pm1,P1)))));

            const vector4double a1p1 = vec_sub(F_1, a2p1);
            const vector4double Kp1  = vec_mul(vec_mul(vec_mul(a1p1,a2p1), vec_sub(vec_mul(G1,vec_add(pp1,P1)),vec_mul(G2,vec_add(pp1,P2)))), myreciprocal<preclevel>(vec_add(vec_mul(vec_mul(a1p1,G2),vec_add(pp1,P2)), vec_mul(vec_mul(a2p1,G1),vec_add(pp1,P1)))));

            vector4double rhs1 = vec_add(Kp1, Km1);

            // column 2
            const vector4double a1m2 = vec_sub(F_1, a2m2);
            const vector4double Km2  = vec_mul(vec_mul(vec_mul(a1m2,a2m2), vec_sub(vec_mul(G1,vec_add(pm2,P1)),vec_mul(G2,vec_add(pm2,P2)))), myreciprocal<preclevel>(vec_add(vec_mul(vec_mul(a1m2,G2),vec_add(pm2,P2)), vec_mul(vec_mul(a2m2,G1),vec_add(pm2,P1)))));

            const vector4double a1p2 = vec_sub(F_1, a2p2);
            const vector4double Kp2  = vec_mul(vec_mul(vec_mul(a1p2,a2p2), vec_sub(vec_mul(G1,vec_add(pp2,P1)),vec_mul(G2,vec_add(pp2,P2)))), myreciprocal<preclevel>(vec_add(vec_mul(vec_mul(a1p2,G2),vec_add(pp2,P2)), vec_mul(vec_mul(a2p2,G1),vec_add(pp2,P1)))));

            vector4double rhs2 = vec_add(Kp2, Km2);

            // column 3
            const vector4double a1m3 = vec_sub(F_1, a2m3);
            const vector4double Km3  = vec_mul(vec_mul(vec_mul(a1m3,a2m3), vec_sub(vec_mul(G1,vec_add(pm3,P1)),vec_mul(G2,vec_add(pm3,P2)))), myreciprocal<preclevel>(vec_add(vec_mul(vec_mul(a1m3,G2),vec_add(pm3,P2)), vec_mul(vec_mul(a2m3,G1),vec_add(pm3,P1)))));

            const vector4double a1p3 = vec_sub(F_1, a2p3);
            const vector4double Kp3  = vec_mul(vec_mul(vec_mul(a1p3,a2p3), vec_sub(vec_mul(G1,vec_add(pp3,P1)),vec_mul(G2,vec_add(pp3,P2)))), myreciprocal<preclevel>(vec_add(vec_mul(vec_mul(a1p3,G2),vec_add(pp3,P2)), vec_mul(vec_mul(a2p3,G1),vec_add(pp3,P1)))));

            vector4double rhs3 = vec_add(Kp3, Km3);
#else /* Does check for semi-positivity of \alpha_k and \rho_k\c_k^2 */
            // column 0
            const vector4double Km0  = _computeK_lambda(a2m0, pm0, G1, G2, P1, P2);
            const vector4double Kp0  = _computeK_lambda(a2p0, pp0, G1, G2, P1, P2);

            vector4double rhs0 = vec_add(Kp0, Km0);

            // column 1
            const vector4double Km1  = _computeK_lambda(a2m1, pm1, G1, G2, P1, P2);
            const vector4double Kp1  = _computeK_lambda(a2p1, pp1, G1, G2, P1, P2);

            vector4double rhs1 = vec_add(Kp1, Km1);

            // column 2
            const vector4double Km2  = _computeK_lambda(a2m2, pm2, G1, G2, P1, P2);
            const vector4double Kp2  = _computeK_lambda(a2p2, pp2, G1, G2, P1, P2);

            vector4double rhs2 = vec_add(Kp2, Km2);

            // column 3
            const vector4double Km3  = _computeK_lambda(a2m3, pm3, G1, G2, P1, P2);
            const vector4double Kp3  = _computeK_lambda(a2p3, pp3, G1, G2, P1, P2);

            vector4double rhs3 = vec_add(Kp3, Km3);
#endif /* SEMIPOSITIVE_K */

            _DIEGO_TRANSPOSE4(rhs0, rhs1, rhs2, rhs3);

            rhs0 = vec_add(rhs0, vec_lda(0L, output + offset_out));
            rhs1 = vec_add(rhs1, vec_lda(0L, output + offset_out + DP));
            rhs2 = vec_add(rhs2, vec_lda(0L, output + offset_out + 2 * DP));
            rhs3 = vec_add(rhs3, vec_lda(0L, output + offset_out + 3 * DP));

            vec_sta(rhs0, 0L, output + offset_out);
            vec_sta(rhs1, 0L, output + offset_out + DP);
            vec_sta(rhs2, 0L, output + offset_out + 2 * DP);
            vec_sta(rhs3, 0L, output + offset_out + 3 * DP);
        }
}

inline void _qpx_zextraterm_K(Real * const a2m, Real * const a2p, Real * const pm, Real * const pp,  Real * const output,
        const Real g1, const Real g2, const Real pc1, const Real pc2)
{
    const vector4double F_1 = vec_splats(1.0f);
    const vector4double G1  = vec_splats(g1);
    const vector4double G2  = vec_splats(g2);
    const vector4double P1  = vec_splats(pc1);
    const vector4double P2  = vec_splats(pc2);

    for(int iy=0; iy<OutputSOA::NY; ++iy)
    {
        Real * const a2mptr = a2m + iy * TempSOA::PITCH;
        Real * const a2pptr = a2p + iy * TempSOA::PITCH;
        Real * const pmptr = pm + iy * TempSOA::PITCH;
        Real * const ppptr = pp + iy * TempSOA::PITCH;
        Real * const outputptr = output + iy * OutputSOA::PITCH;

        for(int ix=0; ix<OutputSOA::NX; ix += 4)
        {
            const vector4double a2m0 = vec_lda(0L, a2mptr + ix);
            const vector4double pm0  = vec_lda(0L, pmptr + ix);

            const vector4double a2p0 = vec_lda(0L, a2pptr + ix);
            const vector4double pp0  = vec_lda(0L, ppptr + ix);

#ifndef SEMIPOSITIVE_K /* Does NOT check for semi-positivity of \alpha_k and \rho_k\c_k^2 */
            const vector4double a1m0 = vec_sub(F_1, a2m0);
            const vector4double Km0  = vec_mul(vec_mul(vec_mul(a1m0,a2m0), vec_sub(vec_mul(G1,vec_add(pm0,P1)),vec_mul(G2,vec_add(pm0,P2)))), myreciprocal<preclevel>(vec_add(vec_mul(vec_mul(a1m0,G2),vec_add(pm0,P2)), vec_mul(vec_mul(a2m0,G1),vec_add(pm0,P1)))));

            const vector4double a1p0 = vec_sub(F_1, a2p0);
            const vector4double Kp0  = vec_mul(vec_mul(vec_mul(a1p0,a2p0), vec_sub(vec_mul(G1,vec_add(pp0,P1)),vec_mul(G2,vec_add(pp0,P2)))), myreciprocal<preclevel>(vec_add(vec_mul(vec_mul(a1p0,G2),vec_add(pp0,P2)), vec_mul(vec_mul(a2p0,G1),vec_add(pp0,P1)))));
#else /* Does check for semi-positivity of \alpha_k and \rho_k\c_k^2 */
            const vector4double Km0  = _computeK_lambda(a2m0, pm0, G1, G2, P1, P2);
            const vector4double Kp0  = _computeK_lambda(a2p0, pp0, G1, G2, P1, P2);
#endif /* SEMIPOSITIVE_K */

            const vector4double oldvalue = vec_lda(0L, outputptr + ix);

            vec_sta(vec_add(oldvalue, vec_add(Kp0, Km0)), 0L, outputptr + ix);
        }
    }
}

void DivSOA2D_QPX::xrhs(const TempSOA& flux, OutputSOA& rhs) const
{
	_qpx_xrhsadd(const_cast<Real*>(flux.ptr(0,0)), &rhs.ref(0,0));
}

void DivSOA2D_QPX::yrhs(const TempSOA& flux, OutputSOA& rhs) const
{
	_qpx_yrhsadd(const_cast<Real*>(flux.ptr(0,0)), &rhs.ref(0,0));
}

void DivSOA2D_QPX::zrhs(const TempSOA& fback, const TempSOA& fforward, OutputSOA& rhs) const
{
	_qpx_zrhsadd(const_cast<Real*>(fback.ptr(0,0)), const_cast<Real*>(fforward.ptr(0,0)), &rhs.ref(0,0));
}

void DivSOA2D_HLLC_QPX_5eq::xextraterm(const TempSOA& um, const TempSOA& up,
                                       const TempSOA& A2m, const TempSOA& A2p,
                                       const TempSOA& pm, const TempSOA& pp,
                                       const TempSOA& am, const TempSOA& ap, const TempSOA& as,
                                       OutputSOA& divu, OutputSOA& sumA2, OutputSOA& sumK)
{
    // consistent HLLC velocity at interface
    _qpx_xextraterm(const_cast<Real *>(um.ptr(0,0)), const_cast<Real *>(up.ptr(0,0)),
                    const_cast<Real *>(am.ptr(0,0)), const_cast<Real *>(ap.ptr(0,0)), const_cast<Real *>(as.ptr(0,0)), &divu.ref(0,0));

    // average alpha 2 term in cell center based on reconstructed value at cell
    // face
    _qpx_xextraterm(const_cast<Real *>(A2m.ptr(0,0)), const_cast<Real *>(A2p.ptr(0,0)), &sumA2.ref(0,0));
#ifndef _NOK_
    // K*div(u) term
    _qpx_xextraterm_K(const_cast<Real *>(A2m.ptr(0,0)), const_cast<Real *>(A2p.ptr(0,0)),
                      const_cast<Real *>(pm.ptr(0,0)), const_cast<Real *>(pp.ptr(0,0)), &sumK.ref(0,0),
                      m_g1, m_g2, m_pc1, m_pc2);
#endif /* _NOK_ */
}

void DivSOA2D_HLLC_QPX_5eq::yextraterm(const TempSOA& vm, const TempSOA& vp,
                                       const TempSOA& A2m, const TempSOA& A2p,
                                       const TempSOA& pm, const TempSOA& pp,
                                       const TempSOA& am, const TempSOA& ap, const TempSOA& as,
                                       OutputSOA& divu, OutputSOA& sumA2, OutputSOA& sumK)
{
    // consistent HLLC velocity at interface
    _qpx_yextraterm(const_cast<Real *>(vm.ptr(0,0)), const_cast<Real *>(vp.ptr(0,0)),
                    const_cast<Real *>(am.ptr(0,0)), const_cast<Real *>(ap.ptr(0,0)), const_cast<Real *>(as.ptr(0,0)), &divu.ref(0,0));

    // average alpha 2 term in cell center based on reconstructed value at cell
    // face
    _qpx_yextraterm(const_cast<Real *>(A2m.ptr(0,0)), const_cast<Real *>(A2p.ptr(0,0)), &sumA2.ref(0,0));
#ifndef _NOK_
    // K*div(u) term
    _qpx_yextraterm_K(const_cast<Real *>(A2m.ptr(0,0)), const_cast<Real *>(A2p.ptr(0,0)),
                      const_cast<Real *>(pm.ptr(0,0)), const_cast<Real *>(pp.ptr(0,0)), &sumK.ref(0,0),
                      m_g1, m_g2, m_pc1, m_pc2);
#endif /* _NOK_ */
}

void DivSOA2D_HLLC_QPX_5eq::zextraterm(const TempSOA& wm0, const TempSOA& wp0, const TempSOA& wm1, const TempSOA& wp1,
                                       const TempSOA& A2m, const TempSOA& A2p, const TempSOA& pm, const TempSOA& pp,
                                       const TempSOA& am0, const TempSOA& ap0, const TempSOA& am1, const TempSOA& ap1,
                                       const TempSOA& as0, const TempSOA& as1,
                                       OutputSOA& divu, OutputSOA& sumA2, OutputSOA& sumK)
{
    // consistent HLLC velocity at interface
    _qpx_zextraterm(const_cast<Real *>(wm0.ptr(0,0)), const_cast<Real *>(wp0.ptr(0,0)),
                    const_cast<Real *>(wm1.ptr(0,0)), const_cast<Real *>(wp1.ptr(0,0)),
                    const_cast<Real *>(am0.ptr(0,0)), const_cast<Real *>(ap0.ptr(0,0)),
                    const_cast<Real *>(am1.ptr(0,0)), const_cast<Real *>(ap1.ptr(0,0)),
                    const_cast<Real *>(as0.ptr(0,0)), const_cast<Real *>(as1.ptr(0,0)),
                    &divu.ref(0,0));

    // average alpha 2 term in cell center based on reconstructed value at cell
    // face
    _qpx_zextraterm(const_cast<Real *>(A2m.ptr(0,0)), const_cast<Real *>(A2p.ptr(0,0)), &sumA2.ref(0,0));
#ifndef _NOK_
    // K*div(u) term
    _qpx_zextraterm_K(const_cast<Real *>(A2m.ptr(0,0)), const_cast<Real *>(A2p.ptr(0,0)),
                      const_cast<Real *>(pm.ptr(0,0)), const_cast<Real *>(pp.ptr(0,0)), &sumK.ref(0,0),
                      m_g1, m_g2, m_pc1, m_pc2);
#endif /* _NOK_ */
}
