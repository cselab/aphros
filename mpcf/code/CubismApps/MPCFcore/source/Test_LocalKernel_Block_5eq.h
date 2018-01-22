/*
 *  Test_LocalKernel_Block_5eq.h
 *  MPCFcore
 *
 *  Created by Ursula Rasthofer on 2016/03/01.
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */

#ifndef TEST_LOCALKERNEL_BLOCK_5EQ_H
#define TEST_LOCALKERNEL_BLOCK_5EQ_H

#include "Test_Utils_5eq.h"

#include "Update.h"

class Test_LocalKernel_Block_5eq
{
  ArgumentParser& parser;

  public:
  
  Test_LocalKernel_Block_5eq(ArgumentParser& par):parser(par){}

  template<typename TKERNEL>
  void accuracy(TKERNEL& kernel)
  {
     // initialization
     Block * blockgold = new Block;
     Block * block = new Block;
     // get reference data from file
     parser.set_strict_mode();
     string fname = parser("-reference-data").asString();
     parser.unset_strict_mode();
     const bool genref = parser("-gen-ref-data").asBool(false);
     if(!genref)
       _initialize_block(*blockgold, fname, true);
     else
       _initialize_block(*blockgold); // get some simple values, test will fail anyways in this case
     _initialize_block(*block, parser, false);

     //run kernel to test
     kernel.compute(&(*block)(0,0,0).dsdt.a1r1, &(*block)(0,0,0).s.a1r1, sizeof(GP)/sizeof(Real));
     if(genref)
       write_block_to_file(*block,"reference_data_update_kernel.txt", true);

     // compute errors
     double tolabs[7]; 
     tolabs[0] = parser("-tolabs-a1r1").asDouble(1.0e-10);
     tolabs[1] = parser("-tolabs-a2r2").asDouble(1.0e-10);
     tolabs[2] = parser("-tolabs-ru").asDouble(1.0e-10);
     tolabs[3] = parser("-tolabs-rv").asDouble(1.0e-10);
     tolabs[4] = parser("-tolabs-rw").asDouble(1.0e-10);
     tolabs[5] = parser("-tolabs-E").asDouble(1.0e-10);
     tolabs[6] = parser("-tolabs-a2").asDouble(1.0e-10);
     const double tolrel = parser("-tolrel").asDouble(1.0e-10);     

     if (!genref)
     {
       Real * const data= &(*block)(0,0,0).s.a1r1;
       Real * const gold_data= &(*blockgold)(0,0,0).s.a1r1;

       const int srcfloats  = sizeof(GP)/sizeof(Real);
       const bool verbose =false;
       if(!verbose)
         printAccuracyTitle();
       for (int i = 0; i < _BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*srcfloats; i += srcfloats)
         check_err(tolabs, tolrel, &data[i], &gold_data[i], 7, verbose);
      }

     delete block;
     delete blockgold;

    return;
  }

  template<typename TKERNEL, typename TREF>
  void accuracy(TKERNEL& kernel, TREF& refkernel)
  {
     // initialization
     Block * blockgold = new Block;
     Block * block = new Block;
     _initialize_block(*blockgold, parser, false);
     _initialize_block(*block, parser, false);

     // run reference kernel
     refkernel.compute(&(*blockgold)(0,0,0).dsdt.a1r1, &(*blockgold)(0,0,0).s.a1r1, sizeof(GP)/sizeof(Real));
     //run kernel to test
     kernel.compute(&(*block)(0,0,0).dsdt.a1r1, &(*block)(0,0,0).s.a1r1, sizeof(GP)/sizeof(Real));

     // compute errors
     double tolabs[7];
     tolabs[0] = parser("-tolabs-a1r1").asDouble(1.0e-10);
     tolabs[1] = parser("-tolabs-a2r2").asDouble(1.0e-10);
     tolabs[2] = parser("-tolabs-ru").asDouble(1.0e-10);
     tolabs[3] = parser("-tolabs-rv").asDouble(1.0e-10);
     tolabs[4] = parser("-tolabs-rw").asDouble(1.0e-10);
     tolabs[5] = parser("-tolabs-E").asDouble(1.0e-10);
     tolabs[6] = parser("-tolabs-a2").asDouble(1.0e-10);
     const double tolrel = parser("-tolrel").asDouble(1.0e-10);
   
     Real * const data= &(*block)(0,0,0).s.a1r1;
     Real * const gold_data= &(*blockgold)(0,0,0).s.a1r1;

     const int srcfloats  = sizeof(GP)/sizeof(Real);
     const bool verbose =false;
     if(!verbose)
       printAccuracyTitle();
     for (int i = 0; i < _BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*srcfloats; i += srcfloats)
       check_err(tolabs, tolrel, &data[i], &gold_data[i], 7, verbose);

     delete block;
     delete blockgold;

    return;
  }



//TODO: performance and profiling not changed in March 2016

                template<typename TKernel> double _benchmarkUP(TKernel kernel, const int NBLOCKS, const int NTIMES)
                {
                        Block * block = new Block[NBLOCKS];

                        for(int i=0; i<NBLOCKS; i++)
                                _initialize_block(block[i]);

                        Timer timer;

                        if (NBLOCKS == 1)
                                kernel.compute(&block[0](0,0,0).dsdt.a1r1, &block[0](0,0,0).s.a1r1, sizeof(GP)/sizeof(Real));

                        timer.start();
                        for(int i=0; i<NTIMES; i++)
                                kernel.compute(&block[i%NBLOCKS](0,0,0).dsdt.a1r1, &block[i%NBLOCKS](0,0,0).s.a1r1, sizeof(GP)/sizeof(Real));

                        const double tCOMPUTE = timer.stop();

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

					for(int i=0; i<NBLOCKS; i++)
						_initialize_block(block[i]);

					Timer timer;

					//run gold
					if (NBLOCKS==1)
						refkernel.compute(&blockgold[0](0,0,0).dsdt.a1r1, &blockgold[0](0,0,0).s.a1r1, sizeof(GP)/sizeof(Real));

#pragma omp barrier
					timer.start();
					for(int i=0; i<NTIMES; i++)
						refkernel.compute(&blockgold[i%NBLOCKS](0,0,0).dsdt.a1r1, &blockgold[i%NBLOCKS](0,0,0).s.a1r1, sizeof(GP)/sizeof(Real));
					const double tgold = timer.stop();

#pragma omp barrier
					const double t = _benchmarkUP(kernel, NBLOCKS, NTIMES);

					delete [] block;
					delete [] blockgold;

#pragma omp critical
					{
						tGOLD += tgold;
						tCOMPUTE += t;
					}
				}

				tCOMPUTE /= COUNT;
				tGOLD /= COUNT;

				kernel.printflops(PEAKPERF * 1e9, PEAKBAND * 1e9, 1, 1, NTIMES, tCOMPUTE, false);

				printf("\tGAIN-OVER-GOLD: %.2fX\n", tGOLD/tCOMPUTE);
			}


		template<typename TKERNEL>
			void profile(TKERNEL& kernel, const double PEAKPERF, const double PEAKBAND, const int NBLOCKS=8*8*8, const int NTIMES=100, bool bAwk=false)
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
					//Block * block = new Block[NBLOCKS];

					//for(int i=0; i<NBLOCKS; i++)
					//	_initialize(block[i]);

					Timer timer;

#pragma omp barrier
					const double t = _benchmarkUP(kernel, NBLOCKS, NTIMES);

					//delete [] block;

#pragma omp critical
					{
						tCOMPUTE += t;
					}
				}

				tCOMPUTE /= COUNT;

				kernel.printflops(PEAKPERF * 1e9, PEAKBAND  * 1e9, 1, 1, NTIMES, tCOMPUTE, false);
			}

};

#endif /* TEST_LOCALKERNEL_BLOCK_5EQ_H */
