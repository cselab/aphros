/*
 *  DivTensor_CPP.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 2/24/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */
#pragma once
#include <vector>
#include <cstdio>
#include <cmath>

#include "SOA2D.h"

typedef SOA2D<-1, _BLOCKSIZE_+1, -1, _BLOCKSIZE_+1, float> InputSOAf_ST;
typedef SOA2D<-1, _BLOCKSIZE_+1, -1, _BLOCKSIZE_+1> InputSOA_ST;
typedef RingSOA2D<-1, _BLOCKSIZE_+1, -1, _BLOCKSIZE_+1, 2> RingInputSOA_ST;

typedef SOA2D<0,_BLOCKSIZE_+1, 0,_BLOCKSIZE_+1, float> TempSOAf_ST;
typedef SOA2D<0,_BLOCKSIZE_+1, 0,_BLOCKSIZE_+1> TempSOA_ST;
typedef RingSOA2D<0, _BLOCKSIZE_+1, 0,_BLOCKSIZE_+1, 2> RingTempSOA_ST;

typedef SOA2D<0,_BLOCKSIZE_+1, 0,_BLOCKSIZE_, float> TempPiXSOAf_ST;
typedef SOA2D<0,_BLOCKSIZE_, 0,_BLOCKSIZE_+1, float> TempPiYSOAf_ST;
typedef SOA2D<0,_BLOCKSIZE_, 0,_BLOCKSIZE_, float> TempPiZSOAf_ST;
typedef SOA2D<0,_BLOCKSIZE_+1, 0,_BLOCKSIZE_> TempPiXSOA_ST;
typedef SOA2D<0,_BLOCKSIZE_, 0,_BLOCKSIZE_+1> TempPiYSOA_ST;
typedef SOA2D<0,_BLOCKSIZE_, 0,_BLOCKSIZE_> TempPiZSOA_ST;
typedef RingSOA2D<0, _BLOCKSIZE_, 0,_BLOCKSIZE_, 2> RingTempPiZSOA_ST;

class DivTensor_CPP
{
public:

	const Real a; //factor that affects the rhs, related to low-storage rk3 (see _copyback method)
	const Real sigma; //amplification factor of the rhs
	const Real h, dtinvh; //grid spacing and "lambda"

protected:

    struct AssumedType { Real r1, r2, u, v, w, E, A2, dummy; };

    RingInputSOA_ST ringu, ringv, ringw, ringls; //slices for the primitive values
    RingTempSOA_ST ringnx, ringny, ringnz; //slices for the corner-gradients

	TempPiXSOA_ST txx, txy, txz, utx; //slices for the tensors in the x-faces
	TempPiYSOA_ST tyx, tyy, tyz, uty; //slices for the tensors in the y-faces
	RingTempPiZSOA_ST ringtzx, ringtzy, ringtzz, ringutz; //slices for the tensors in the z-faces

	OutputSOA rhsu, rhsv, rhsw, rhss; //output slices

    virtual void _convert(const Real * const gptfirst, const int gptfloats, const int rowgpts);

    virtual void _corners(const InputSOA_ST& ls0, const InputSOA_ST& ls1, TempSOA_ST& nx, TempSOA_ST& ny, TempSOA_ST& nz);

	virtual void _tensor_xface(const TempSOA_ST& nx0, const TempSOA_ST& nx1,
							   const TempSOA_ST& ny0, const TempSOA_ST& ny1,
							   const TempSOA_ST& nz0, const TempSOA_ST& nz1);

	virtual void _tensor_yface(const TempSOA_ST& nx0, const TempSOA_ST& nx1,
							   const TempSOA_ST& ny0, const TempSOA_ST& ny1,
							   const TempSOA_ST& nz0, const TempSOA_ST& nz1);

	virtual void _tensor_zface(const TempSOA_ST& nx, const TempSOA_ST& ny, const TempSOA_ST& nz,
							   TempPiZSOA_ST& tzx, TempPiZSOA_ST& tzy, TempPiZSOA_ST& tzz);

	virtual void _udot_tx(const InputSOA_ST& u, const InputSOA_ST& v, const InputSOA_ST& w);
	virtual void _udot_ty(const InputSOA_ST& u, const InputSOA_ST& v, const InputSOA_ST& w);

	virtual void _udot_tz(const InputSOA_ST& u0, const InputSOA_ST& v0, const InputSOA_ST& w0,
						  const InputSOA_ST& u1, const InputSOA_ST& v1, const InputSOA_ST& w1,
						  const TempPiZSOA_ST& tzx, const TempPiZSOA_ST& tzy, const TempPiZSOA_ST& tzz,
						  TempPiZSOA_ST& utz);

	virtual void _div_dxy();
	virtual void _div_dz(const TempPiZSOA_ST& tzx0, const TempPiZSOA_ST& tzy0, const TempPiZSOA_ST& tzz0, const TempPiZSOA_ST& utz0,
						 const TempPiZSOA_ST& tzx1, const TempPiZSOA_ST& tzy1, const TempPiZSOA_ST& tzz1, const TempPiZSOA_ST& utz1);

    virtual void _copyback(Real * const gptfirst, const int gptfloats, const int rowgpts);

	void _input_next() { ringu.next(); ringv.next(); ringw.next(); ringls.next(); }

    void _normals_next() { ringnx.next(); ringny.next(); ringnz.next(); }

    void _tensors_next() { ringtzx.next(); ringtzy.next(); ringtzz.next(); ringutz.next(); }

public:

	DivTensor_CPP(const Real a = 1, const Real dtinvh = 1, const Real h = 1, const Real sigma=1):
    a(a), dtinvh(dtinvh), h(h), sigma(sigma)
	{
	}

	//main loop, processing the z-slices in sequence. all the logic is in here
    void compute(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
                 Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts)
	{
		_convert(srcfirst, srcfloats, rowsrcs);
		_input_next();

		_convert(srcfirst + srcfloats*slicesrcs, srcfloats, rowsrcs);

		_corners(ringls(-1), ringls(0), ringnx.ref(), ringny.ref(), ringnz.ref());
		_tensor_zface(ringnx(), ringny(), ringnz(), ringtzx.ref(), ringtzy.ref(), ringtzz.ref());
		_udot_tz(ringu(-1), ringv(-1), ringw(-1),  ringu(0), ringv(0), ringw(0), ringtzx(), ringtzy(), ringtzz(), ringutz.ref());

		for(int islice=0; islice<_BLOCKSIZE_; islice++)
		{
			_tensors_next();
			_input_next();
			_normals_next();

			_convert(srcfirst + (islice+2)*srcfloats*slicesrcs, srcfloats, rowsrcs);
			_corners(ringls(-1), ringls(0), ringnx.ref(), ringny.ref(), ringnz.ref());

			_tensor_xface(ringnx(0), ringnx(1), ringny(0), ringny(1), ringnz(0), ringnz(1));
			_tensor_yface(ringnx(0), ringnx(1), ringny(0), ringny(1), ringnz(0), ringnz(1));
			_tensor_zface(ringnx(), ringny(), ringnz(), ringtzx.ref(), ringtzy.ref(), ringtzz.ref());

			_udot_tx(ringu(-1), ringv(-1), ringw(-1));
			_udot_ty(ringu(-1), ringv(-1), ringw(-1));
			_udot_tz(ringu(-1), ringv(-1), ringw(-1), ringu(0), ringv(0), ringw(0), ringtzx(), ringtzy(), ringtzz(), ringutz.ref());

			_div_dxy();
			_div_dz(ringtzx(-1), ringtzy(-1), ringtzz(-1), ringutz(-1), ringtzx(0), ringtzy(0), ringtzz(0), ringutz(0));

			_copyback(dstfirst + islice*dstfloats*slicedsts, dstfloats, rowdsts);
		}
	}

	virtual void hpc_info(float& flop_convert, int& traffic_convert,
						  float& flop_corners, int& traffic_corner,
						  float& flop_tensorface, int& traffic_tensorface,
						  float& flop_tdotu, int& traffic_tdotu,
						  float& flop_div, int& traffic_div,
						  float& flop_copyback, int& traffic_copyback,
						  size_t& footprint)
	{
		const int ninputs = (int)powf(_BLOCKSIZE_ + 2, 3);
		const int ncorners = (int)powf(_BLOCKSIZE_ + 1, 3);
		const int ntensors = 3 * (_BLOCKSIZE_ + 1) * _BLOCKSIZE_ * _BLOCKSIZE_;
		const int ncells = (int)powf(_BLOCKSIZE_, 3);

		flop_convert = 3 * ninputs;
		traffic_convert = (4 + 4) * sizeof(Real) * ninputs;
		flop_corners = 21 * ncorners;
		traffic_corner = (8 + 3) * sizeof(Real) * ncorners;
		flop_tensorface = 23 * ntensors;
		traffic_tensorface = (12 + 3) * sizeof(Real) * ntensors;
		flop_tdotu = 8 * ntensors;
		traffic_tdotu = (9 + 1) * sizeof(Real) * ntensors;
		flop_div = 5 * 4 * ncells;
		traffic_div = (6 * 4 + 4) * sizeof(Real) * ncells;
		flop_copyback = 3 * 4 * ncells;
		traffic_copyback = (4 + 4) * sizeof(Real) * ncells;

		footprint = sizeof(*this);
	}

	void printflops(const float PEAKPERF_CORE, const float PEAKBAND, const int NCORES, const int NT, const int NBLOCKS, const float MEASUREDTIME_, const bool bAwk=false)
	{
		const float MEASUREDTIME = MEASUREDTIME_;
		const float PEAKPERF = PEAKPERF_CORE*NCORES;

		float flop_convert, flop_corners, flop_tensorface, flop_tdotu, flop_div, flop_copyback;
		int traffic_convert, traffic_corner, traffic_tensorface, traffic_tdotu, traffic_div, traffic_copyback;
		size_t footprint;

		hpc_info(flop_convert, traffic_convert, flop_corners, traffic_corner, flop_tensorface, traffic_tensorface,
				 flop_tdotu, traffic_tdotu, flop_div, traffic_div, flop_copyback, traffic_copyback, footprint);

		double texpected_ai = 0;
		double totflop = 0;

		//compute texpected_ai
		{
			std::vector<float> ai(6), flop(6);

			ai[0] = flop_convert/traffic_convert;
			ai[1] = flop_corners/traffic_corner;
			ai[2] = flop_tensorface/traffic_tensorface;
			ai[3] = flop_tdotu/traffic_tdotu;
			ai[4] = flop_div/traffic_div;
			ai[5] = flop_copyback/traffic_copyback;

			flop[0] = flop_convert;
			flop[1] = flop_corners;
			flop[2] = flop_tensorface;
			flop[3] = flop_tdotu;
			flop[4] = flop_div;
			flop[5] = flop_copyback;

			for(int i=0; i<ai.size(); ++i)
				texpected_ai += NT * NBLOCKS * flop[i] / min(PEAKPERF, PEAKBAND*ai[i]);

			for(int i=0; i<ai.size(); ++i)
				totflop += NT * NBLOCKS * flop[i];
		}

		const double ai_overall = min((double)PEAKPERF , (totflop / texpected_ai) / PEAKBAND);

		const double inout_footprint =  NT * NBLOCKS * (5 * sizeof(Real) * powf(_BLOCKSIZE_+2, 3) + 2 * 4 * sizeof(Real) * powf(_BLOCKSIZE_, 3));

		const double oi_overall = totflop/(inout_footprint + NT * (2 + 1) * footprint);
		const double texpected_oi = totflop/min((double)PEAKPERF, PEAKBAND*oi_overall);

		const double perf_measured = 1e-9*totflop/MEASUREDTIME;

		printPerformanceTitle();
		printf("\tINTERMEDIATE MEMORY FOOTPRINT: %.4f MB\tTOTAL TRAFFIC: %.4f MB\n", footprint/1024./1024, (NT * NBLOCKS * footprint + inout_footprint)/1024./1024);
		printf("\tASSUMING PP: %.2f GFLOP/s (PER CORE), %.2f GFLOP/s (OVERALL)\n\tPB: %.2f GB/s (OVERALL)\n", PEAKPERF_CORE*1e-9, PEAKPERF*1e-9, PEAKBAND*1e-9);
		printf("\tRIDGE AT %.2f FLOP/B\n", PEAKPERF/PEAKBAND);
		printf("\tDIVTENSOR THIS ONE IS %.2f GFLOP/s,\t\"per block\" %.2f FLOP/B [AI] - %.2f FLOP/B [OI]\n", perf_measured, ai_overall, oi_overall);
		printf("\tTIME PER BLOCK: %.5f ms (expected %.5f [AI] - %.5f [OI] ms)\n",  1e3*MEASUREDTIME/(NT * NBLOCKS), 1e3*texpected_ai/(NT * NBLOCKS), 1e3*texpected_oi/(NT * NBLOCKS));
		printf("\tExpected Performance is: %.2f  GFLOP/s [AI], %.2f  GFLOP/s [OI]\n", totflop*1e-9/texpected_ai, totflop*1e-9/texpected_oi);
		printf("\tEFFICIENCY: %.2f%% [AI] - %.2f%% [OI], HW-UTILIZATION: %.2f%%\n", 100.*min(1., texpected_ai/MEASUREDTIME), 100.*texpected_oi/MEASUREDTIME, 100*perf_measured*1e9/PEAKPERF);
		printEndLine();
	}
};
