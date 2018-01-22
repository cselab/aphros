#pragma once

#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include <algorithm>

#include "FFTWBaseMPI.h"

using namespace std;

class FFTWScalarMPI : public FFTWBaseMPI
{

public:

  FFTWScalarMPI(const int desired_threads): FFTWBaseMPI(desired_threads), initialized(false) {  }

  virtual void dispose()
  {
    if (initialized)
    {
      initialized = false;

#ifndef _FFTW_SP_COMP_
      fftw_destroy_plan(this->fwd);

      fftw_destroy_plan(bwd);

      fftw_free(data);
#else
      fftwf_destroy_plan(fwd);
      fftwf_destroy_plan(bwd);

      fftwf_free(data);
#endif
      FFTWBaseMPI::dispose();
    }

    return;
  }

protected:

  virtual void _setup(MYReal *& rhs , const size_t nx, const size_t ny, const size_t nz)
  {
    if (!initialized)
    {
      initialized = true;
#ifndef _FFTW_SP_COMP_
      alloc_local = fftw_mpi_local_size_3d_transposed(nx, ny, nz/2+1, MPI::COMM_WORLD, &local_n0, &local_0_start, &local_n1, &local_1_start);

      rhs = fftw_alloc_real(2*alloc_local);

      fwd = fftw_mpi_plan_dft_r2c_3d(nx, ny, nz, rhs, (mycomplex *)rhs, MPI::COMM_WORLD, FFTW_MPI_TRANSPOSED_OUT | FFTW_MEASURE);
      bwd = fftw_mpi_plan_dft_c2r_3d(nx, ny, nz, (mycomplex *)rhs, rhs, MPI::COMM_WORLD, FFTW_MPI_TRANSPOSED_IN | FFTW_MEASURE);
#else
      alloc_local = fftwf_mpi_local_size_3d_transposed(nx, ny, nz/2+1, MPI::COMM_WORLD, &local_n0, &local_0_start, &local_n1, &local_1_start);

      rhs = fftwf_alloc_real(2*alloc_local);

      fwd = fftwf_mpi_plan_dft_r2c_3d(nx, ny, nz, rhs, (mycomplex *)rhs, MPI::COMM_WORLD, FFTW_MPI_TRANSPOSED_OUT | FFTW_MEASURE);
      bwd = fftwf_mpi_plan_dft_c2r_3d(nx, ny, nz, (mycomplex *)rhs, rhs, MPI::COMM_WORLD, FFTW_MPI_TRANSPOSED_IN | FFTW_MEASURE);
#endif
    }
    return;
  }

  virtual void fft_backward()
  {
    if (initialized)
    {
#ifndef _FFTW_SP_COMP_
      fftw_execute(bwd);
#else
      fftwf_execute(bwd);
#endif
    }
    else
    {
      cout << "FFTW scalar not initialized " << endl;
      abort();
    }
    return;
  }

  virtual void fft_forward()
  {
    if (initialized)
    {
#ifndef _FFTW_SP_COMP_
      fftw_execute(fwd);
#else
      fftwf_execute(fwd);
#endif
    }
    else
    {
      cout << "FFTW scalar not initialized " << endl;
      abort();
    }
    return;
  }

protected:

  bool initialized;

  myplan fwd, bwd;
  ptrdiff_t alloc_local, local_n0, local_0_start, local_n1, local_1_start;
  MYReal * data; // rhs in _setup, out in cub2fftw and fftw2cub

};
