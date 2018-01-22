/*
 *  MaxSpeedOfSound.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 6/15/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <cstdio>

#include "common.h"
#include "SOA2D.h"

class MaxSpeedOfSound_CPP
{
public:
	virtual Real compute(const Real * const src, const int gptfloats) const;

	static void printflops(const float PEAKPERF_CORE, const float PEAKBAND, const int NCORES, const int NT, const int NBLOCKS, float MEASUREDTIME, const bool bAwk=false)
	{
		const float PEAKPERF = PEAKPERF_CORE*NCORES;

		//FLOP estimation
		const double flopcount = 5 + 4 + mysqrt_flops<preclevel>() + 10 + 3 * myminmax_flops() + 2 * myreciprocal_flops<preclevel>();
		const double GFLOPSOS = NBLOCKS*_BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*flopcount*1.e-9;

		//FLOP/s estimation
		const float OISOS = flopcount/(8 * sizeof(Real));
		const float AISOS = flopcount/(8 * sizeof(Real));
		const double EPERFSOS = min(OISOS*PEAKBAND, PEAKPERF);
		const double EPERFSOSAI = min(AISOS*PEAKBAND, PEAKPERF);

		const double TSOS = 1.e9*GFLOPSOS/EPERFSOS;

		printPerformanceTitle();
		printf("\tSOS TOTAL GFLOPS: %.2f\n", GFLOPSOS);
		printf("\tSOS ASSUMING PP: %.2f GFLOP/s (PER CORE), %.2f GFLOP/s (OVERALL)\n\tPB: %.2f GB/s (OVERALL)\n", PEAKPERF_CORE*1e-9, PEAKPERF*1e-9, PEAKBAND*1e-9);
		printf("\tPER ITERATION: %d FLOPs\n", (int)flopcount);
		printf("\tSOS RIDGE AT %.2f FLOP/B\n", PEAKPERF/PEAKBAND);
		printf("\tSOS THIS ONE IS %.2f GFLOP/s,\t\"PER SOS\" %.2f FLOP/B\n", GFLOPSOS/MEASUREDTIME, OISOS);
		printf("\tSOS TIME PER BLOCK: %.5f ms (expected %.5f ms)\n",  1e3*MEASUREDTIME/NBLOCKS, 1e3*TSOS/NBLOCKS);
		printf("\tSOS Expected Performance is: %.2f GFLOP/s [AI], %.2f GFLOP/s [OI]\n", GFLOPSOS/(1.e9*GFLOPSOS/EPERFSOSAI), GFLOPSOS/TSOS);
		printf("\tSOS EFFICIENCY: %.2f%% [AI] - %.2f%% [OI], HW-UTILIZATION: %.2f%%\n", 100.*(1.e9*GFLOPSOS/EPERFSOSAI)/MEASUREDTIME, 100.*TSOS/MEASUREDTIME, 100*(GFLOPSOS/MEASUREDTIME*1e9)/PEAKPERF);
		printEndLine();
	}
};
