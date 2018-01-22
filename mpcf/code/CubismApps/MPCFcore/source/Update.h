/*
 *  Update.h
 *  MPCFcore
 *
 *  Created by Babak Hejazialhosseini  on 6/9/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once
#include <cassert>
#include <cstdlib>
#include <cstdio>

#include "common.h"

class Update_CPP
{
protected:

	Real m_b;

	inline bool _is_aligned(const void * const ptr, unsigned int alignment) const
	{
		return ((size_t)ptr) % alignment == 0;
	}

public:

	Update_CPP(Real b=1): m_b(b) {}

	void compute(const Real * const src, Real * const dst, const int gptfloats) const;

	static void printflops(const float PEAKPERF_CORE, const float PEAKBAND, const size_t NCORES, const size_t NT, const size_t NBLOCKS, const float MEASUREDTIME, const bool bAwk=false)
	{
		const float PEAKPERF = PEAKPERF_CORE*NCORES;

		//FLOP estimation
		const double GFLOPUPDATE = 8 * 2.e-9 * NBLOCKS * (_BLOCKSIZE_ * _BLOCKSIZE_ * _BLOCKSIZE_);

		//FLOP/s estimation
		const float OIUpdate = 2.f/(3 * sizeof(Real));
		const float AIUpdate = 2.f/(4 * sizeof(Real));
		const double EPERFUPDATE   = min(OIUpdate*   PEAKBAND, PEAKPERF);
		const double EPERFUPDATEAI   = min(AIUpdate*   PEAKBAND, PEAKPERF);

		//execution time estimation
		const double TUPDATE = 1.e9*GFLOPUPDATE/EPERFUPDATE;
		const double TUPDATEAI = 1.e9*GFLOPUPDATE/EPERFUPDATEAI;

		printPerformanceTitle();
		printf("\tUP TOTAL GFLOPS: %.2f\n", GFLOPUPDATE);
		printf("\tUP ASSUMING PP: %.2f GFLOP/s (PER CORE), %.2f GFLOP/s (OVERALL)\n\tPB: %.2f GB/s (OVERALL)\n", PEAKPERF_CORE*1e-9, PEAKPERF*1e-9, PEAKBAND*1e-9);
		printf("\tUP RIDGE AT %.2f FLOP/B\n", PEAKPERF/PEAKBAND);
		printf("\tUP THIS ONE IS %.2f GFLOP/s,\t\"PER Update\" %.2f FLOP/B\n", GFLOPUPDATE/MEASUREDTIME, OIUpdate);
		printf("\tUP TIME PER BLOCK: %.5f ms (expected %.5f ms)\n",  1e3*MEASUREDTIME/NBLOCKS, 1e3*TUPDATE/NBLOCKS);
		printf("\tUP Expected Performance is: %.2f GFLOP/s [AI], %.2f GFLOP/s [OI]\n", GFLOPUPDATE/TUPDATEAI, GFLOPUPDATE/TUPDATE);
		printf("\tUP EFFICIENCY: %.2f%% [AI] - %.2f%% [OI], HW-UTILIZATION: %.2f%%\n", 100.*TUPDATEAI/MEASUREDTIME, 100.*TUPDATE/MEASUREDTIME, 100*(GFLOPUPDATE/MEASUREDTIME*1e9)/PEAKPERF);
		printEndLine();
	}
};
