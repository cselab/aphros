/*
 *  Test_LocalKernel_Val_5eq.h
 *  MPCFcore
 *
 *  Created by Ursula Rasthofer on 2016/03/01.
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */

#ifndef TEST_LOCALKERNEL_VAL_5EQ_H
#define TEST_LOCALKERNEL_VAL_5EQ_H

#include "Test_Utils_5eq.h"
#include "check_errors.h"
#include <Timer.h>

class Test_LocalKernel_Val_5eq
{

  ArgumentParser& parser;

  public:

  Test_LocalKernel_Val_5eq(ArgumentParser& par):parser(par){}


  // testing againts reference values
  template<typename TKERNEL>
  void accuracy(TKERNEL& kernel)
  {
    // initialization
    Block * block = new Block;
   _initialize_block(*block, parser, true);

    // compute max sos for this block
    // run kernel
    const Real result = kernel.compute(&(*block)(0,0,0).s.a1r1, sizeof(GP)/sizeof(Real));

    // get reference value
    parser.set_strict_mode();
    const Real refval = parser("-reference-value").asDouble();
    parser.unset_strict_mode();

    // check results
    const double tolabs = parser("-tolabs").asDouble(1.0e-10);
    const double tolrel = parser("-tolrel").asDouble(1.0e-10);
    check_err(tolabs, tolrel, refval, result, parser("-verb").asBool(true));

    delete block;

    return;
  }


  // comparing two kernels
  template<typename TTESTKERNEL, typename TREFKERNEL>
  void accuracy(TTESTKERNEL& kernel, TREFKERNEL& refkernel)
  {
    // initialization
    Block * block = new Block;
    _initialize_block(*block, parser, true);

    // compute max sos for this block
    // run reference kernel
    const Real refval = refkernel.compute(&(*block)(0,0,0).s.a1r1, sizeof(GP)/sizeof(Real));
    // run kernel to test
    const Real result = kernel.compute(&(*block)(0,0,0).s.a1r1, sizeof(GP)/sizeof(Real));

    // check result
    const double tolabs = parser("-tolabs").asDouble(1.0e-10);
    const double tolrel = parser("-tolrel").asDouble(1.0e-10);
    check_err(tolabs, tolrel, refval, result);

    delete block;

    return;
  }



//TODO: performance and profiling not changed in March 2016

		template<typename TKernel> double _benchmarkSOS(TKernel kernel, const int NBLOCKS, const int NTIMES)
		{
			Block * block = new Block[NBLOCKS];

			for(int i=0; i<NBLOCKS; i++)
				_initialize_block(block[i]);

			Timer timer;

			if (NBLOCKS == 1)
				kernel.compute(&block[0](0,0,0).s.a1r1, sizeof(GP)/sizeof(Real));

#pragma omp master
			{
				HPM_Start("SOSDIEGO");
			}
#pragma omp barrier
			timer.start();
			for(int i=0; i<NTIMES; i++)
				kernel.compute(&block[i%NBLOCKS](0,0,0).s.a1r1, sizeof(GP)/sizeof(Real));

			const double tCOMPUTE = timer.stop();
#pragma omp barrier
#pragma omp master
			{
				HPM_Stop("SOSDIEGO");
			}
			delete [] block;

			return tCOMPUTE;
		}

		template<typename TTESTKERNEL, typename TREFKERNEL>
			void performance(TTESTKERNEL& kernel, TREFKERNEL& refkernel, const double PEAKPERF, const double PEAKBAND, const int NBLOCKS=8*8*8, const int NTIMES=100, bool bAwk=false)
			{
				int COUNT = 0;

#pragma omp parallel
				{
#pragma omp critical
					{
						COUNT++;
					}
				}

				//measure performance
				float tGOLD = 0, tCOMPUTE = 0;

#pragma omp parallel
				{

					Block * blockgold = new Block[NBLOCKS];
					Block * block = new Block[NBLOCKS];

					for(int i=0; i<NBLOCKS; i++)
						_initialize_block(blockgold[i]);

					if (NBLOCKS == 1)
						refkernel.compute(&blockgold[0](0,0,0).s.a1r1, sizeof(GP)/sizeof(Real));

					//run gold
#pragma omp barrier
					Timer timer;

					timer.start();
					for(int i=0; i<NTIMES; i++)
						refkernel.compute(&blockgold[i%NBLOCKS](0,0,0).s.a1r1, sizeof(GP)/sizeof(Real));
					const double tgold = timer.stop();

#pragma omp barrier
					const double t = _benchmarkSOS(kernel, NBLOCKS, NTIMES);

#pragma omp critical
					{
						tGOLD += tgold;
						tCOMPUTE += t;
					}

					delete [] block;
					delete [] blockgold;
				}

				tCOMPUTE /= COUNT;
				tGOLD /= COUNT;

				kernel.printflops(PEAKPERF * 1e9, PEAKBAND * 1e9, 1, 1, NTIMES, tCOMPUTE, false);

				printf("\tGAIN-OVER-GOLD: %.2fX\n", tGOLD/tCOMPUTE);
			}


		template<typename TKERNEL>
			void profile(TKERNEL& kernel, const double PEAKPERF = 2.66*8/(sizeof(Real)/4)*1e9, const double PEAKBAND = 4.5*1e9, const int NBLOCKS=8*8*8, const int NTIMES=100, bool bAwk=false)
			{
				int COUNT = 0;

#pragma omp parallel
				{
#pragma omp critical
					{
						COUNT++;
					}
				}

				//measure performance
				float tCOMPUTE = 0;

#pragma omp parallel
				{
					//	Block * block = new Block[NBLOCKS];

					//	for(int i=0; i<NBLOCKS; i++)
					//		_initialize(block[i]);

					Timer timer;

#pragma omp barrier
					const double t = _benchmarkSOS(kernel, NBLOCKS, NTIMES);

					//	delete [] block;

#pragma omp critical
					{
						tCOMPUTE += t;
					}
				}

				tCOMPUTE /= COUNT;

				kernel.printflops(PEAKPERF * 1e9, PEAKBAND * 1e9, 1, 1, NTIMES, tCOMPUTE, false);
			}
};

#endif /* TEST_LOCALKERNEL_VAL_5EQ_H */
