/*
 *  Test_LocalKernel.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 6/17/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "MaxSpeedOfSound.h"
#include "Update.h"

class Test_LocalKernel
{

	protected:

		Real gamma1, gamma2, smoothlength;

		virtual void _initialize(Block& block)
		{
			//printf("initilaizigin from random data....\n");
			srand48(61651);

			for(int iz = 0; iz<_BLOCKSIZE_; iz++)
				for(int iy = 0; iy<_BLOCKSIZE_; iy++)
					for(int ix = 0; ix<_BLOCKSIZE_; ix++)
					{
						const int a = iy + 3;
						const int b = iz + 3;
						const int c = ix + 3;
						const double L = _BLOCKSIZE_;

						block(ix, iy, iz).clear();
						block(ix, iy, iz).s.r = (40+a*a+(c+3)*L + (a+3)*L*L)/(double)L ;//drand48()+iz;
						block(ix, iy, iz).s.u =(40+a*a+(c+3)*L + (a+3)*L*L)/(double)L/(double)L;// drand48()+ix;
						block(ix, iy, iz).s.v = (40+a*a+(c+3)*L + (a+3)*L*L)/(double)L/(double)L;//drand48()+iy;
						block(ix, iy, iz).s.w = (40+a*a+(c+3)*L + (a+3)*L*L)/(double)L/(double)L;//drand48()+(ix*iy/(double) _BLOCKSIZE_);
						block(ix, iy, iz).s.s =  100+b/(double)L;///drand48()+1+(iy*iz)/(double) _BLOCKSIZE_;
						block(ix, iy, iz).s.G = (1+c+b*L + a*L*L)/(double)L;
						block(ix, iy, iz).s.P = (1+c+b*L + a*L*L)/(double)L;//block(ix, iy, iz).s.levelset = -1 +  (iz+i)x+iy)/(double) _BLOCKSIZE_;

						block(ix, iy, iz).dsdt.r = drand48()+iz;
						block(ix, iy, iz).dsdt.u = drand48()+ix;
						block(ix, iy, iz).dsdt.v = drand48()+iy;
						block(ix, iy, iz).dsdt.w = drand48()+(ix*iy/(double) _BLOCKSIZE_);
						block(ix, iy, iz).dsdt.s = drand48()+1+(iy*iz)/(double) _BLOCKSIZE_;
						//block(ix, iy, iz).dsdt.levelset = -1 +  (iz+ix+iy)/(double) _BLOCKSIZE_;
					}
		}

		virtual void _initialize_from_file(Block& block)
		{
			FILE * f = fopen("test_rank0", "r");

			int c = 0;
			while (!feof(f))
			{
				float r,u,v,w,s,l;
				const int values = fscanf(f, "%e %e %e %e %e %e\n", &r, &u, &v, &w, &s, &l);
				assert(values==6);

				const int ix = (c % 16);
				const int iy = (c/16 % 16) ;
				const int iz = (c/16/16);
				block(ix, iy, iz).clear();
				block(ix, iy, iz).s.r = r;
				block(ix, iy, iz).s.u = u;
				block(ix, iy, iz).s.v = v;
				block(ix, iy, iz).s.w = w;
				//block(ix, iy, iz).s.levelset = s;
				block(ix, iy, iz).s.s = l;

				c++;
			}

			fclose(f);
		}

	public:

		template<typename TSOS>
			void accuracy(TSOS& kernel, MaxSpeedOfSound_CPP& refkernel, double accuracy=1e-4, bool bAwk=false)
			{
				Block * blockgold = new Block;
				Block * block = new Block;

				_initialize(*blockgold);
				_initialize(*block);

				const Real v1 = refkernel.compute(&(*blockgold)(0,0,0).s.r, sizeof(GP)/sizeof(Real));

				//run Kernel
				const Real v2 = kernel.compute(&(*blockgold)(0,0,0).s.r, sizeof(GP)/sizeof(Real));

				printAccuracyTitle();
				printf("\tVALUES: %e (reference) %e (tested kernel)\n", v1, v2);

				printf("\tERROR: %e (relative: %e)\n", v1-v2, (v1-v2)/max(fabs(v1), max(fabs(v2), (Real)accuracy)));

				delete block;
				delete blockgold;
			}

		template<typename TUPDATE>
			void accuracy(TUPDATE& kernel, Update_CPP& refkernel, double accuracy=1e-4, bool bAwk=false)
			{
				Block * blockgold = new Block;
				Block * block = new Block;

				_initialize(*blockgold);
				_initialize(*block);

				refkernel.compute(&(*blockgold)(0,0,0).dsdt.r, &(*blockgold)(0,0,0).s.r, sizeof(GP)/sizeof(Real));

				//run Kernel
				kernel.compute(&(*block)(0,0,0).dsdt.r, &(*block)(0,0,0).s.r, sizeof(GP)/sizeof(Real));

				printAccuracyTitle();
				//block->compare(*blockgold, accuracy, typeid(TUPDATE).name(), false);
				Real * const data= &(*block)(0,0,0).s.r;
				Real * const gold_data= &(*blockgold)(0,0,0).s.r;

				const int srcfloats  = sizeof(GP)/sizeof(Real);
				for (int i = 0; i < _BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*srcfloats; i += srcfloats)
				{
					check_error(accuracy, &data[i], &gold_data[i], 7);
				}

				delete block;
				delete blockgold;
			}

		template<typename TKernel> double _benchmarkSOS(TKernel kernel, const int NBLOCKS, const int NTIMES)
		{
			Block * block = new Block[NBLOCKS];

			for(int i=0; i<NBLOCKS; i++)
				_initialize(block[i]);

			Timer timer;

			if (NBLOCKS == 1)
				kernel.compute(&block[0](0,0,0).s.r, sizeof(GP)/sizeof(Real));

#pragma omp master
			{
				HPM_Start("SOSDIEGO");
			}
#pragma omp barrier
			timer.start();
			for(int i=0; i<NTIMES; i++)
				kernel.compute(&block[i%NBLOCKS](0,0,0).s.r, sizeof(GP)/sizeof(Real));

			const double tCOMPUTE = timer.stop();
#pragma omp barrier
#pragma omp master
			{
				HPM_Stop("SOSDIEGO");
			}
			delete [] block;

			return tCOMPUTE;
		}

		template<typename TKernel> double _benchmarkUP(TKernel kernel, const int NBLOCKS, const int NTIMES)
		{
			Block * block = new Block[NBLOCKS];

			for(int i=0; i<NBLOCKS; i++)
				_initialize(block[i]);

			Timer timer;

			if (NBLOCKS == 1)
				kernel.compute(&block[0](0,0,0).dsdt.r, &block[0](0,0,0).s.r, sizeof(GP)/sizeof(Real));

			timer.start();
			for(int i=0; i<NTIMES; i++)
				kernel.compute(&block[i%NBLOCKS](0,0,0).dsdt.r, &block[i%NBLOCKS](0,0,0).s.r, sizeof(GP)/sizeof(Real));

			const double tCOMPUTE = timer.stop();

			delete [] block;

			return tCOMPUTE;
		}


		template<typename TSOS>
			void performance(TSOS& kernel, MaxSpeedOfSound_CPP& refkernel, const double PEAKPERF, const double PEAKBAND, const int NBLOCKS=8*8*8, const int NTIMES=100, bool bAwk=false)
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
						_initialize(blockgold[i]);

					if (NBLOCKS == 1)
						refkernel.compute(&blockgold[0](0,0,0).s.r, sizeof(GP)/sizeof(Real));

					//run gold
#pragma omp barrier
					Timer timer;

					timer.start();
					for(int i=0; i<NTIMES; i++)
						refkernel.compute(&blockgold[i%NBLOCKS](0,0,0).s.r, sizeof(GP)/sizeof(Real));
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

		template<typename TUPDATE>
			void performance(TUPDATE& kernel, Update_CPP& refkernel, const double PEAKPERF, const double PEAKBAND, const int NBLOCKS=8*8*8, const int NTIMES=100, bool bAwk=false)
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
						_initialize(blockgold[i]);

					for(int i=0; i<NBLOCKS; i++)
						_initialize(block[i]);

					Timer timer;

					//run gold
					if (NBLOCKS==1)
						refkernel.compute(&blockgold[0](0,0,0).dsdt.r, &blockgold[0](0,0,0).s.r, sizeof(GP)/sizeof(Real));

#pragma omp barrier
					timer.start();
					for(int i=0; i<NTIMES; i++)
						refkernel.compute(&blockgold[i%NBLOCKS](0,0,0).dsdt.r, &blockgold[i%NBLOCKS](0,0,0).s.r, sizeof(GP)/sizeof(Real));
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


		template<typename TUPDATE>
			void profile_update(TUPDATE& kernel, const double PEAKPERF, const double PEAKBAND, const int NBLOCKS=8*8*8, const int NTIMES=100, bool bAwk=false)
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

		template<typename TSOS>
			void profile_maxsos(TSOS& kernel, const double PEAKPERF = 2.66*8/(sizeof(Real)/4)*1e9, const double PEAKBAND = 4.5*1e9, const int NBLOCKS=8*8*8, const int NTIMES=100, bool bAwk=false)
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
