/*
 *  FFTWBaseMPI.h
 *
 *  Created by Ursula Rasthofer on 6/29/15.
 *  Copyright 2015 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <fftw3.h>
#ifndef _FFTW_SP_COMP_
typedef fftw_complex mycomplex;
typedef fftw_plan myplan;
typedef double MYReal;
#else
typedef fftwf_complex mycomplex;
typedef fftwf_plan myplan;
typedef float MYReal;
#endif
#include "fftw3-mpi.h"

#include <iostream>
#include <stdlib.h>
using namespace std;

class FFTWBaseMPI
{
  static int registered_objects_;
  static bool initialized_; //a la singleton

//protected:
//  myplan fwd, bwd;
//  ptrdiff_t alloc_local, local_n0, local_0_start, local_n1, local_1_start;
//  MYReal * data; // rhs in _setup, out in cub2fftw and fftw2cub
//private:
  static void _setup(const int desired_threads)
  {
      if (!initialized_)
      {
          initialized_ = true;

          const int supported_threads = MPI::Query_thread();

          if (supported_threads>=MPI_THREAD_FUNNELED)
          {
#ifndef _FFTW_SP_COMP_
              const int retval = fftw_init_threads();
#else
              const int retval = fftwf_init_threads();
#endif
              if(retval==0)
              {
                  cout << "FFTWBaseMPI::setup(): Oops the call to fftw_init_threads() returned zero. Aborting\n";
                  abort();
              }
              else
              {
#ifndef _FFTW_SP_COMP_
                fftw_plan_with_nthreads(desired_threads);
#else
                fftwf_plan_with_nthreads(desired_threads);
#endif
              }
          }

#ifndef _FFTW_SP_COMP_
          fftw_mpi_init();
#else
          fftwf_mpi_init();
#endif
      }

      registered_objects_++;
  }

public:

  FFTWBaseMPI(const int desired_threads) { _setup(desired_threads); }

  static void dispose()
  {
      registered_objects_--;

      if (registered_objects_ == 0)
      {
#ifndef _FFTW_SP_COMP_
          fftw_mpi_cleanup();
#else
          fftwf_mpi_cleanup();
#endif
      }
  }
};
