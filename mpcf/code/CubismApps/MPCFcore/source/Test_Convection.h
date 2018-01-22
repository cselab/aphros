/*
 *  Convection_Test.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 5/17/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <typeinfo>

#include <Timer.h>

#include "common.h"
#include "TestTypes.h"
#include "Convection_CPP.h"
#include "check_errors.h"
using namespace std;

class Test_Convection
{
	virtual void _initialize_lab(TestLab& lab);
	virtual void _initialize_block(Block& block);
	void _print(Block& block);

	Real dtinvh;

public:

	inline void _gold(Convection_CPP &kernel, TestLab& lab, Block& block)
	{
		const Real * const srcfirst = &(lab)(-3,-3, -3).s.r;
		const int srcfloats  = sizeof(GP)/sizeof(Real);
		const int rowsrcs = _BLOCKSIZE_+6;
		const int slicesrcs =  (_BLOCKSIZE_+6)*(_BLOCKSIZE_+6);

		Real * const dstfirst = &(block)(0,0,0).dsdt.r;
		const int dstfloats =  sizeof(GP)/sizeof(Real);
		const int rowdsts =  _BLOCKSIZE_;
		const int slicedsts = _BLOCKSIZE_*_BLOCKSIZE_;

		kernel.compute(srcfirst, srcfloats, rowsrcs, slicesrcs, dstfirst, dstfloats, rowdsts, slicedsts);
	}


	template<typename TCOV>
	inline void _apply_kernel(TCOV &kernel, TestLab& lab, Block& block)
	{
		const Real * const srcfirst = &(lab)(-3,-3, -3).s.r;
		const int srcfloats  = sizeof(GP)/sizeof(Real);
		const int rowsrcs = _BLOCKSIZE_+6;
		const int slicesrcs =  (_BLOCKSIZE_+6)*(_BLOCKSIZE_+6);

		Real * const dstfirst = &(block)(0,0,0).dsdt.r;
		const int dstfloats =  sizeof(GP)/sizeof(Real);
		const int rowdsts =  _BLOCKSIZE_;
		const int slicedsts = _BLOCKSIZE_*_BLOCKSIZE_;

		kernel.compute(srcfirst, srcfloats, rowsrcs, slicesrcs, dstfirst, dstfloats, rowdsts, slicedsts);
	}


	template<typename TCOV>
	void accuracy(TCOV& kernel, double accuracy=1e-4, bool bAwk=false)
	{
		Convection_CPP refkernel(0,1);
		dtinvh = kernel.dtinvh;

		TestLab * lab = new TestLab;
		Block * blockgold = new Block;
		Block * block = new Block;

		_initialize_lab(*lab);
		_initialize_block(*blockgold);
		_initialize_block(*block);

		// run Reference Kernel (Convection_CPP)
		#if 0
		{
			const Real * const srcfirst = &(*lab)(-3,-3, -3).s.r;
			const int srcfloats  = sizeof(GP)/sizeof(Real);
			const int rowsrcs = _BLOCKSIZE_+6;
			const int slicesrcs =  (_BLOCKSIZE_+6)*(_BLOCKSIZE_+6);

			Real * const dstfirst = &(*blockgold)(0,0,0).dsdt.r;
			const int dstfloats =  sizeof(GP)/sizeof(Real);
			const int rowdsts =  _BLOCKSIZE_;
			const int slicedsts = _BLOCKSIZE_*_BLOCKSIZE_;

			refkernel.compute(srcfirst, srcfloats, rowsrcs, slicesrcs, dstfirst, dstfloats, rowdsts, slicedsts);

		}
		#else
			//_gold(refkernel, *lab, *blockgold);
			/* _apply_kernel(refkernel, *lab, *blockgold); */
		#endif

                //run Kernel
		#if 0
		{
			const Real * const srcfirst = &(*lab)(-3,-3, -3).s.r;
			const int srcfloats  = sizeof(GP)/sizeof(Real);
			const int rowsrcs = _BLOCKSIZE_+6;
			const int slicesrcs =  (_BLOCKSIZE_+6)*(_BLOCKSIZE_+6);

			Real * const dstfirst = &(*block)(0,0,0).dsdt.r;
			const int dstfloats =  sizeof(GP)/sizeof(Real);
			const int rowdsts =  _BLOCKSIZE_;
			const int slicedsts = _BLOCKSIZE_*_BLOCKSIZE_;

			kernel.compute(srcfirst, srcfloats, rowsrcs, slicesrcs, dstfirst, dstfloats, rowdsts, slicedsts);

		}
		#else
			_apply_kernel(kernel, *lab, *block);
		#endif

		printAccuracyTitle();
		string kernelname = typeid(kernel).name();

		// check errors
		{
			/* Real * const data= &(*block)(0,0,0).dsdt.r; */
			/* Real * const gold_data= &(*blockgold)(0,0,0).dsdt.r; */

			/* const int srcfloats  = sizeof(GP)/sizeof(Real); */
			/* for (int i = 0; i < _BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*srcfloats; i += srcfloats) */
			/* { */
			/* 	check_error(accuracy, &data[i], &gold_data[i], 7); */
			/* } */

			//check_error(accuracy, data, gold_data,_BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*srcfloats); // sizeof(GP)/sizeof(Real));

		}

		printEndLine();

		delete block;
		delete blockgold;
		delete lab;
        }

protected:
	template<typename FS> double _benchmark(FS fs, const int NBLOCKS, const int NTIMES)
	{
		TestLab * lab = new TestLab[NBLOCKS];
		Block * block = new Block[NBLOCKS];

		//printf("_benchmark: NBLOCKS = %d NTIMES = %d\n", NBLOCKS, NTIMES);

		for(int i=0; i<NBLOCKS; i++) {
			_initialize_lab(lab[i]);
			_initialize_block(block[i]);
		}

		double tCOMPUTE = 0;

		//measure performance
		{
			Timer timer;

			timer.start();

			//run FS
			const int srcfloats  = sizeof(GP)/sizeof(Real);
			const int rowsrcs = _BLOCKSIZE_+6;
			const int slicesrcs =  (_BLOCKSIZE_+6)*(_BLOCKSIZE_+6);

			const int dstfloats =  sizeof(GP)/sizeof(Real);
			const int rowdsts =  _BLOCKSIZE_;
			const int slicedsts = _BLOCKSIZE_*_BLOCKSIZE_;

			if (NBLOCKS == 1) {
				const Real * const srcfirst = &(lab[0])(-3,-3, -3).s.r;
				Real * const dstfirst = &(block[0])(0,0,0).dsdt.r;

				fs.compute(srcfirst, srcfloats, rowsrcs, slicesrcs, dstfirst, dstfloats, rowdsts, slicedsts);
			}

			#pragma omp barrier
			timer.start();
			for(int i=0; i<NTIMES; i++) {
				const Real * const srcfirst = &(lab[i%NBLOCKS])(-3,-3, -3).s.r;
				Real * const dstfirst = &(block[i%NBLOCKS])(0,0,0).dsdt.r;
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

	template<typename FS> void performance(FS& fs, const double PEAKPERF = 2.66*8/(sizeof(Real)/4)*1e9, const double PEAKBAND = 4.5*1e9, const int NBLOCKS=8*8*8, const int NTIMES=100, bool bAwk=false)
	{
		Convection_CPP refkernel(0,1);
		dtinvh = fs.dtinvh;

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
			TestLab * lab = new TestLab[NBLOCKS];
			Block * blockgold = new Block[NBLOCKS];

			for(int i=0; i<NBLOCKS; i++) {
				_initialize_lab(lab[i]);
				_initialize_block(blockgold[i]);
			}

			Timer timer;

			//run gold
			if (NBLOCKS == 1)
				_gold(refkernel, lab[0], blockgold[0]);

#pragma omp barrier

			timer.start();
			for(int i=0; i<NTIMES; i++)
				_gold(refkernel, lab[i%NBLOCKS], blockgold[i%NBLOCKS]);
			const double tgold = timer.stop();

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
		dtinvh = fs.dtinvh;

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
		dtinvh = fs.dtinvh;

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

