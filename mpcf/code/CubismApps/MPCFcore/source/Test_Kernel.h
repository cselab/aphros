/*
 *  Test_Kernel.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 5/17/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#ifndef TEST_KERNEL_H
#define TEST_KERNEL_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <typeinfo>

#include <Timer.h>
#include "common.h"
#include "Test_Utils_5eq.h"
#include "TestTypes.h"
#include "check_errors.h"

using namespace std;

template<int WIDTH> // stencil width changes with kernels
class Test_Kernel
{
  ArgumentParser& parser;

public:

  Test_Kernel(ArgumentParser& par) : parser(par) {return;};

  template<typename TKernel>
  inline void _apply_kernel(TKernel &kernel, TestLab<WIDTH>& lab, Block& block)
  {
    const Real * const srcfirst = &(lab)(-WIDTH,-WIDTH, -WIDTH).s.a1r1;
    const int srcfloats  = sizeof(GP)/sizeof(Real);
    const int rowsrcs = _BLOCKSIZE_+2*WIDTH;
    const int slicesrcs =  (_BLOCKSIZE_+2*WIDTH)*(_BLOCKSIZE_+2*WIDTH);

    Real * const dstfirst = &(block)(0,0,0).dsdt.a1r1;
    const int dstfloats =  sizeof(GP)/sizeof(Real);
    const int rowdsts =  _BLOCKSIZE_;
    const int slicedsts = _BLOCKSIZE_*_BLOCKSIZE_;

    kernel.compute(srcfirst, srcfloats, rowsrcs, slicesrcs, dstfirst, dstfloats, rowdsts, slicedsts);
  }

  template<typename TKERNEL>
  void accuracy(TKERNEL& kernel)
  {
     // initialization
     TestLab<WIDTH> * lab = new TestLab<WIDTH>;
     Block * blockgold = new Block;
     Block * block = new Block;
     // get reference data from file
     parser.set_strict_mode();
     string fname = parser("-reference-data").asString();
     parser.unset_strict_mode();
     _initialize_lab<WIDTH>(*lab, parser);
     const bool genref = parser("-gen-ref-data").asBool(false);
     if(!genref)
       _initialize_block(*blockgold, fname, false);
     else
       _initialize_block(*blockgold);
     _initialize_block(*block);

     //run kernel to test
     _apply_kernel(kernel, *lab, *block);
     if(genref)
       write_block_to_file(*block,"reference_data_kernel.txt",false);

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
       Real * const data= &(*block)(0,0,0).dsdt.a1r1;
       Real * const gold_data= &(*blockgold)(0,0,0).dsdt.a1r1;

       const int srcfloats  = sizeof(GP)/sizeof(Real);
       const bool verbose = parser("verb").asBool(true);
       if(verbose)
         printAccuracyTitle();
       for (int i = 0; i < _BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*srcfloats; i += srcfloats)
         check_err(tolabs, tolrel, &data[i], &gold_data[i], 7, verbose);
     }

     delete block;
     delete blockgold;
     delete lab;

     return;
  }

  template<typename TTESTKERNEL, typename TREFKERNEL>
  void accuracy(TTESTKERNEL& kernel, TREFKERNEL& refkernel)
  {
     // initialization
     TestLab<WIDTH> * lab = new TestLab<WIDTH>;
     Block * blockgold = new Block;
     Block * block = new Block;
     _initialize_lab<WIDTH>(*lab,parser);
     _initialize_block(*blockgold);
     _initialize_block(*block);

     // run reference kernel
     _apply_kernel(refkernel, *lab, *blockgold);
     // run kernel to test
     _apply_kernel(kernel, *lab, *block);

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

     Real * const data= &(*block)(0,0,0).dsdt.a1r1;
     Real * const gold_data= &(*blockgold)(0,0,0).dsdt.a1r1;

     const int srcfloats  = sizeof(GP)/sizeof(Real);
     const bool verbose = parser("verb").asBool(true);
     if(verbose)
       printAccuracyTitle();
     for (int i = 0; i < _BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*srcfloats; i += srcfloats)
       check_err(tolabs, tolrel, &data[i], &gold_data[i], 7, verbose);

     delete block;
     delete blockgold;
     delete lab;

     return;
  }



protected:
	template<typename FS> double _benchmark(FS fs, const int NBLOCKS, const int NTIMES)
	{
		TestLab<WIDTH> * lab = new TestLab<WIDTH>[NBLOCKS];
		Block * block = new Block[NBLOCKS];

		//printf("_benchmark: NBLOCKS = %d NTIMES = %d\n", NBLOCKS, NTIMES);

		for(int i=0; i<NBLOCKS; i++) {
			_initialize_lab<WIDTH>(lab[i]);
			_initialize_block(block[i]);
		}

		double tCOMPUTE = 0;

		//measure performance
		{
			Timer timer;

			timer.start();

			//run FS
			const int srcfloats  = sizeof(GP)/sizeof(Real);
			const int rowsrcs = _BLOCKSIZE_+2*WIDTH;
			const int slicesrcs =  (_BLOCKSIZE_+2*WIDTH)*(_BLOCKSIZE_+2*WIDTH);

			const int dstfloats =  sizeof(GP)/sizeof(Real);
			const int rowdsts =  _BLOCKSIZE_;
			const int slicedsts = _BLOCKSIZE_*_BLOCKSIZE_;

			if (NBLOCKS == 1) {
				const Real * const srcfirst = &(lab[0])(-WIDTH,-WIDTH, -WIDTH).s.a1r1;
				Real * const dstfirst = &(block[0])(0,0,0).dsdt.a1r1;

				fs.compute(srcfirst, srcfloats, rowsrcs, slicesrcs, dstfirst, dstfloats, rowdsts, slicedsts);
			}

			#pragma omp barrier
			timer.start();
			for(int i=0; i<NTIMES; i++) {
				const Real * const srcfirst = &(lab[i%NBLOCKS])(-WIDTH,-WIDTH, -WIDTH).s.a1r1;
				Real * const dstfirst = &(block[i%NBLOCKS])(0,0,0).dsdt.a1r1;
				fs.compute(srcfirst, srcfloats, rowsrcs, slicesrcs, dstfirst, dstfloats, rowdsts, slicedsts);
			}
			tCOMPUTE = timer.stop();
			//printf("tCOMPUTE = %lf ms\n", tCOMPUTE*1e3);
		}

		delete [] lab;
		delete [] block;

		return tCOMPUTE;
	}

public:

	template<typename FS, typename REFKERNEL> void performance(FS& fs, REFKERNEL refkernel, const double PEAKPERF = 2.66*8/(sizeof(Real)/4)*1e9, const double PEAKBAND = 4.5*1e9, const int NBLOCKS=8*8*8, const int NTIMES=100, bool bAwk=false)
	{
		//printf("CONV: performance\n");
		int COUNT = 0;

#pragma omp parallel
		{
#pragma omp critical
			{
				COUNT++;
			}
		}

		double tGOLD = 0, tCOMPUTE = 0;

#pragma omp parallel
		{
			TestLab<WIDTH> * lab = new TestLab<WIDTH>[NBLOCKS];
			Block * blockgold = new Block[NBLOCKS];

			for(int i=0; i<NBLOCKS; i++) {
				_initialize_lab<WIDTH>(lab[i]);
				_initialize_block(blockgold[i]);
			}

			Timer timer;

			//run gold
			if (NBLOCKS == 1)
				_apply_kernel(refkernel, lab[0], blockgold[0]);

#pragma omp barrier

            // TODO: (fabianw@mavt.ethz.ch; Thu 18 Aug 2016 11:42:14 AM CEST)
            // For identical kernels, the gain is not 1.00X if the code below
            // is used.
			// timer.start();
			// for(int i=0; i<NTIMES; i++)
			// 	_apply_kernel(refkernel, lab[i%NBLOCKS], blockgold[i%NBLOCKS]);
			// const double tgold = timer.stop();
			const double tgold = _benchmark(refkernel, NBLOCKS, NTIMES);

#pragma omp barrier

			const double t = _benchmark(fs, NBLOCKS, NTIMES);

#pragma omp barrier

#pragma omp critical
			{
				tGOLD += tgold;
				tCOMPUTE += t;
			}
		}

		tCOMPUTE /= COUNT;
		tGOLD /= COUNT;


		fs.printflops(PEAKPERF, PEAKBAND, 1, NTIMES, 1, tCOMPUTE);

		printf("\tGAIN-OVER-GOLD: %.2fX\n", tGOLD/tCOMPUTE);
	}

	template<typename FS> void profile(FS& fs, const double PEAKPERF = 2.66*8/(sizeof(Real)/4)*1e9, const double PEAKBAND = 4.5*1e9, const int NBLOCKS=8*8*8, const int NTIMES=100, bool bAwk=false)
	{

		int COUNT = 0;

		//printf(":profile\n");

#pragma omp parallel
		{
#pragma omp critical
			{
				COUNT++;
			}
		}

		double tCOMPUTE = 0;

#pragma omp parallel
		{
//			TestLab * lab = new TestLab[NBLOCKS];
//			Block * block = new Block[NBLOCKS];

//			for(int i=0; i<NBLOCKS; i++) {
//				_initialize_lab(lab[i]);
//				_initialize_block(block[i]);
//			}

			//float tCOMPUTE = 0;

			const double t = _benchmark(fs, NBLOCKS, NTIMES);

#pragma omp critical
			{
				tCOMPUTE += t;
			}

//			delete block;
//			delete lab;
		}

		tCOMPUTE /= COUNT;

		fs.printflops(PEAKPERF, PEAKBAND, 1, NTIMES, 1, tCOMPUTE);
	}

	template<typename FS> void profile2(FS& fs, const double PEAKPERF = 2.66*8/(sizeof(Real)/4)*1e9, const double PEAKBAND = 4.5*1e9, const int NBLOCKS=8*8*8, const int NTIMES=100, bool bAwk=false)
	{

		int COUNT = 0;

		//printf(":profile2\n");

//#pragma omp parallel
		{
#pragma omp critical
			{
				COUNT++;
			}
		}

		double tCOMPUTE = 0;

//#pragma omp parallel
		{
//			TestLab * lab = new TestLab[NBLOCKS];
//			Block * block = new Block[NBLOCKS];

//			for(int i=0; i<NBLOCKS; i++) {
//				_initialize_lab(lab[i]);
//				_initialize_block(block[i]);
//			}

			//float tCOMPUTE = 0;

			const double t = _benchmark(fs, NBLOCKS, NTIMES);

#pragma omp critical
			{
				tCOMPUTE += t;
			}

//			delete block;
//			delete lab;
		}

		tCOMPUTE /= COUNT;

		fs.printflops(PEAKPERF, PEAKBAND, 1, NTIMES, 1, tCOMPUTE);
	}
};

#endif /* TEST_KERNEL_H */
