#pragma once

#include "FFTWScalarMPI.h"

class FFTWVectorMPI : public FFTWScalarMPI
{

protected:

  void _setup(MYReal *& rhs, const size_t nx, const size_t ny, const size_t nz)
  {
    if (!(this->initialized))
    {
      this->initialized = true;

      const ptrdiff_t dim[3] = {static_cast<ptrdiff_t>(nx), static_cast<ptrdiff_t>(ny) ,static_cast<ptrdiff_t>(nz)};

#ifndef _FFTW_SP_COMP_
      this->alloc_local = fftw_mpi_local_size_many_transposed(3, dim, 3, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, MPI::COMM_WORLD, &this->local_n0, &this->local_0_start, &this->local_n1, &this->local_1_start);
      rhs = fftw_alloc_real(2*this->alloc_local);
      this->fwd = fftw_mpi_plan_many_dft_r2c(3, dim, 3, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, rhs, (mycomplex *)rhs, MPI::COMM_WORLD, FFTW_MPI_TRANSPOSED_OUT | FFTW_MEASURE);
      this->bwd = fftw_mpi_plan_many_dft_c2r(3, dim, 3, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, (mycomplex *)rhs, rhs, MPI::COMM_WORLD, FFTW_MPI_TRANSPOSED_IN | FFTW_MEASURE);

#else
      this->alloc_local = fftwf_mpi_local_size_many_transposed(3, dim, 3, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, MPI::COMM_WORLD, &this->local_n0, &this->local_0_start, &this->local_n1, &this->local_1_start);
      rhs = fftw_alloc_real(2*this->alloc_local);
      this->fwd = fftwf_mpi_plan_many_dft_r2c(3, dim, 3, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, rhs, (mycomplex *)rhs, MPI::COMM_WORLD, FFTW_MPI_TRANSPOSED_OUT | FFTW_MEASURE);
      this->bwd = fftwf_mpi_plan_many_dft_c2r(3, dim, 3, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, (mycomplex *)rhs, rhs, MPI::COMM_WORLD, FFTW_MPI_TRANSPOSED_IN | FFTW_MEASURE);
#endif
    }
  }

public:
 FFTWVectorMPI(const int desired_threads) : FFTWScalarMPI(desired_threads)
{
  this->initialized = false;
}


};
