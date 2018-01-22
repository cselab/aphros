#pragma once

#include "PoissonSolverScalarFFTW_MPI.h"

template<typename TGrid, typename TStreamer>
class PoissonSolverVectorFFTW_MPI : public PoissonSolverScalarFFTW_MPI< TGrid, TStreamer >
{

protected:
	void _setup(MYReal *& rhs, const size_t nx, const size_t ny, const size_t nz, MPI_Comm comm)
	{
	  if (!(this->initialized))
		{
			this->initialized = true;

			const ptrdiff_t dim[3] = {static_cast<ptrdiff_t>(nx), static_cast<ptrdiff_t>(ny) ,static_cast<ptrdiff_t>(nz)};

#ifndef _FFTW_SP_COMP_
			this->alloc_local = fftw_mpi_local_size_many_transposed(3, dim, TStreamer::channels, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, comm, &this->local_n0, &this->local_0_start, &this->local_n1, &this->local_1_start);
			rhs = fftw_alloc_real(2*this->alloc_local);
			this->fwd = fftw_mpi_plan_many_dft_r2c(3, dim, TStreamer::channels, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, rhs, (mycomplex *)rhs, comm, FFTW_MPI_TRANSPOSED_OUT | FFTW_MEASURE);
			this->bwd = fftw_mpi_plan_many_dft_c2r(3, dim, TStreamer::channels, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, (mycomplex *)rhs, rhs, comm, FFTW_MPI_TRANSPOSED_IN | FFTW_MEASURE);

#else
			this->alloc_local = fftwf_mpi_local_size_many_transposed(3, dim, TStreamer::channels, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, comm, &this->local_n0, &this->local_0_start, &this->local_n1, &this->local_1_start);
			rhs = fftw_alloc_real(2*this->alloc_local);
			this->fwd = fftwf_mpi_plan_many_dft_r2c(3, dim, TStreamer::channels, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, rhs, (mycomplex *)rhs, comm, FFTW_MPI_TRANSPOSED_OUT | FFTW_MEASURE);
			this->bwd = fftwf_mpi_plan_many_dft_c2r(3, dim, TStreamer::channels, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, (mycomplex *)rhs, rhs, comm, FFTW_MPI_TRANSPOSED_IN | FFTW_MEASURE);
#endif
		}
	}

	void _solve(mycomplex * in_out, const size_t nx, const size_t ny, const size_t nz, const size_t nz_hat, const MYReal norm_factor, const MYReal h)
	{
	  if (TStreamer::channels == 1)
	    {
	      //cout << "PoissonSolverVectorFFTW_MPI::PoissonSolverVectorFFTW_MPI(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be greater than 1. Aborting\n";
	      //abort();
         // do nothing, this also works
	    }

		const MYReal h2 = h*h;
		const MYReal factor = h2*norm_factor;
#pragma omp parallel for
		for(int j = 0; j<this->local_n1; ++j)
			for(int i=0; i<nx; ++i)
				for(int k = 0; k<nz_hat; ++k)
				{
					const int linidx = TStreamer::channels*( (j*nx+i)*nz_hat + k );

					//assert(linidx >=0 && linidx<nx*ny*nz_hat*TStreamer::channels); // linking error with openmp

					const MYReal denom = 32.*(cos(2.*M_PI*i/nx)+cos(2.*M_PI*(this->local_1_start+j)/ny)+cos(2.0*M_PI*k/nz))-2.*(cos(4.*M_PI*i/nx)+cos(4.*M_PI*(this->local_1_start+j)/ny)+cos(4.*M_PI*k/nz))-90.;

					const MYReal inv_denom = (denom==0)? 0.:1./denom;
					const MYReal fatfactor = 12. * inv_denom * factor;

					for(int c=0; c<TStreamer::channels; ++c)
					{
						((mycomplex *)in_out)[linidx+c][0] *= fatfactor;
						((mycomplex *)in_out)[linidx+c][1] *= fatfactor;
					}
				}

        //this is sparta!
        if (this->local_1_start == 0)
            for(int c=0; c<TStreamer::channels; ++c)
                ((mycomplex *)in_out)[c][0] = ((mycomplex *)in_out)[c][1] = 0;
	}

	void _solveSpectral(mycomplex * in_out, const size_t nx, const size_t ny, const size_t nz, const size_t nz_hat, const MYReal norm_factor, const MYReal h)
	{
        if (TStreamer::channels == 1)
	    {
            //cout << "PoissonSolverVectorFFTW_MPI::PoissonSolverVectorFFTW_MPI(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be greater than 1. Aborting\n";
            //abort();
            // do nothing, this also works
	    }

        const MYReal waveFactX = 2.0*M_PI/(nx*h);
        const MYReal waveFactY = 2.0*M_PI/(ny*h);
        const MYReal waveFactZ = 2.0*M_PI/(nz*h);

#pragma omp parallel for
		for(int j = 0; j<this->local_n1; ++j)
			for(int i=0; i<nx; ++i)
				for(int k = 0; k<nz_hat; ++k)
				{
					const int linidx = TStreamer::channels*( (j*nx+i)*nz_hat + k );
					//assert(linidx >=0 && linidx<nx*ny*nz_hat); // linking error with openmp

                    const int kx = (i <= nx/2) ? i : -(nx-i);
                    const int ky = (this->local_1_start+j <= ny/2) ? this->local_1_start+j : -(ny-(this->local_1_start+j));
                    const int kz = (k <= nz/2) ? k : -(nz-k);
                    const MYReal rkx = kx*waveFactX;
                    const MYReal rky = ky*waveFactY;
                    const MYReal rkz = kz*waveFactZ;

                    const MYReal kinv = (kx==0 && ky==0 && kz==0) ? 0.0 : -1.0/(rkx*rkx+rky*rky+rkz*rkz);
					for(int c=0; c<TStreamer::channels; ++c)
					{
						((mycomplex *)in_out)[linidx+c][0] *= kinv*norm_factor;
						((mycomplex *)in_out)[linidx+c][1] *= kinv*norm_factor;
					}
				}

        //this is sparta!
        if (this->local_1_start == 0)
            for(int c=0; c<TStreamer::channels; ++c)
                ((mycomplex *)in_out)[c][0] = ((mycomplex *)in_out)[c][1] = 0;
	}

    void _solveVelocity(mycomplex * in_out, const size_t nx, const size_t ny, const size_t nz, const size_t nz_hat, const MYReal norm_factor, const MYReal h)
	{
        if (TStreamer::channels != 3)
	    {
            cout << "PoissonSolverVectorFFTW_MPI::PoissonSolverVectorFFTW_MPI(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 3. Aborting\n";
            abort();
	    }

        const MYReal waveFactX = 2.0*M_PI/(nx*h);
        const MYReal waveFactY = 2.0*M_PI/(ny*h);
        const MYReal waveFactZ = 2.0*M_PI/(nz*h);
#pragma omp parallel for
		for(int j = 0; j<this->local_n1; ++j)
			for(int i=0; i<nx; ++i)
				for(int k = 0; k<nz_hat; ++k)
				{
					const int linidx = TStreamer::channels*( (j*nx+i)*nz_hat + k );
					//assert(linidx >=0 && linidx<nx*ny*nz_hat*TStreamer::channels); // linking error with openmp

					// wave number
                    const int kx = (i <= nx/2) ? i : -(nx-i);
                    const int ky = (this->local_1_start+j <= ny/2) ? this->local_1_start+j : -(ny-(this->local_1_start+j));
                    const int kz = (k <= nz/2) ? k : -(nz-k);
                    const MYReal rkx = kx*waveFactX;
                    const MYReal rky = ky*waveFactY;
                    const MYReal rkz = kz*waveFactZ;
                    const MYReal kinv = (kx==0 && ky==0 && kz==0) ? 0.0 : 1.0/(rkx*rkx+rky*rky+rkz*rkz);
                    const MYReal wR[3] = {in_out[linidx][0],in_out[linidx+1][0],in_out[linidx+2][0]};
                    const MYReal wC[3] = {in_out[linidx][1],in_out[linidx+1][1],in_out[linidx+2][1]};

                    ((mycomplex *)in_out)[linidx+0][0] = (-rky*wC[2] + rkz*wC[1])*kinv*norm_factor;
                    ((mycomplex *)in_out)[linidx+0][1] = ( rky*wR[2] - rkz*wR[1])*kinv*norm_factor;
                    ((mycomplex *)in_out)[linidx+1][0] = (-rkz*wC[0] + rkx*wC[2])*kinv*norm_factor;
                    ((mycomplex *)in_out)[linidx+1][1] = ( rkz*wR[0] - rkx*wR[2])*kinv*norm_factor;
                    ((mycomplex *)in_out)[linidx+2][0] = (-rkx*wC[1] + rky*wC[0])*kinv*norm_factor;
                    ((mycomplex *)in_out)[linidx+2][1] = ( rkx*wR[1] - rky*wR[0])*kinv*norm_factor;

				}

        //this is sparta!
        if (this->local_1_start == 0)
            for(int c=0; c<TStreamer::channels; ++c)
                ((mycomplex *)in_out)[c][0] = ((mycomplex *)in_out)[c][1] = 0;
	}

    void _solveReproject(mycomplex * in_out, const size_t nx, const size_t ny, const size_t nz, const size_t nz_hat, const MYReal norm_factor, const MYReal h)
	{
        if (TStreamer::channels != 3)
	    {
            cout << "PoissonSolverVectorFFTW_MPI::PoissonSolverVectorFFTW_MPI(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 3. Aborting\n";
            abort();
	    }

        const MYReal waveFactX = 2.0*M_PI/(nx*h);
        const MYReal waveFactY = 2.0*M_PI/(ny*h);
        const MYReal waveFactZ = 2.0*M_PI/(nz*h);
#pragma omp parallel for
		for(int j = 0; j<this->local_n1; ++j)
			for(int i=0; i<nx; ++i)
				for(int k = 0; k<nz_hat; ++k)
				{
					const int linidx = TStreamer::channels*( (j*nx+i)*nz_hat + k );
					//assert(linidx >=0 && linidx<nx*ny*nz_hat*TStreamer::channels); // linking error with openmp

					// wave number
                    const int kx = (i <= nx/2) ? i : -(nx-i);
                    const int ky = (this->local_1_start+j <= ny/2) ? this->local_1_start+j : -(ny-(this->local_1_start+j));
                    const int kz = (k <= nz/2) ? k : -(nz-k);
                    const MYReal rkx = kx*waveFactX;
                    const MYReal rky = ky*waveFactY;
                    const MYReal rkz = kz*waveFactZ;
                    const MYReal kinv = (kx==0 && ky==0 && kz==0) ? 0.0 : 1.0/(rkx*rkx+rky*rky+rkz*rkz);
                    const MYReal wR[3] = {in_out[linidx][0],in_out[linidx+1][0],in_out[linidx+2][0]};
                    const MYReal wC[3] = {in_out[linidx][1],in_out[linidx+1][1],in_out[linidx+2][1]};
                    const MYReal wdotkR = (wR[0]*rkx + wR[1]*rky + wR[2]*rkz)*kinv;
                    const MYReal wdotkC = (wC[0]*rkx + wC[1]*rky + wC[2]*rkz)*kinv;
                    const MYReal rk[3] = {rkx,rky,rkz};
					for(int c=0; c<TStreamer::channels; ++c)
					{
						((mycomplex *)in_out)[linidx+c][0] = (wR[c] - wdotkR*rk[c])*norm_factor;
						((mycomplex *)in_out)[linidx+c][1] = (wC[c] - wdotkC*rk[c])*norm_factor;
					}
				}

        //this is sparta!
        if (this->local_1_start == 0)
            for(int c=0; c<TStreamer::channels; ++c)
                ((mycomplex *)in_out)[c][0] = ((mycomplex *)in_out)[c][1] = 0;
	}
public:
 PoissonSolverVectorFFTW_MPI(const int desired_threads) : PoissonSolverScalarFFTW_MPI< TGrid, TStreamer >(desired_threads)
{
  this->initialized = false;
}

	void solveVelocity(TGrid& grid)
	{
		const int bs[3] = {TGrid::BlockType::sizeX, TGrid::BlockType::sizeY, TGrid::BlockType::sizeZ};
		const size_t gsize[3] = {grid.getBlocksPerDimension(0)*bs[0], grid.getBlocksPerDimension(1)*bs[1], grid.getBlocksPerDimension(2)*bs[2]};

        MPI_Comm comm = grid.getCartComm();

		_setup(this->data, gsize[0], gsize[1], gsize[2], comm);

		this->_cub2fftw(grid, this->data, this->local_n0, gsize[1], gsize[2], gsize[2]/2+1);

#ifndef _FFTW_SP_COMP_
		fftw_execute(this->fwd);
#else
		fftwf_execute(this->fwd);
#endif

		const MYReal norm_factor = 1./(gsize[0]*gsize[1]*gsize[2]);
		const MYReal h = grid.getBlocksInfo().front().h_gridpoint;

        _solveVelocity((mycomplex *)this->data, gsize[0], gsize[1], gsize[2], gsize[2]/2+1, norm_factor, h);

#ifndef _FFTW_SP_COMP_
		fftw_execute(this->bwd);
#else
		fftwf_execute(this->bwd);
#endif

		this->_fftw2cub(this->data, grid, gsize[0], gsize[1], gsize[2], gsize[2]/2+1);

		//profiler.printSummary();
	}

	void reproject(TGrid& grid)
	{
		const int bs[3] = {TGrid::BlockType::sizeX, TGrid::BlockType::sizeY, TGrid::BlockType::sizeZ};
		const size_t gsize[3] = {grid.getBlocksPerDimension(0)*bs[0], grid.getBlocksPerDimension(1)*bs[1], grid.getBlocksPerDimension(2)*bs[2]};

        MPI_Comm comm = grid.getCartComm();

		_setup(this->data, gsize[0], gsize[1], gsize[2], comm);

		this->_cub2fftw(grid, this->data, this->local_n0, gsize[1], gsize[2], gsize[2]/2+1);

#ifndef _FFTW_SP_COMP_
		fftw_execute(this->fwd);
#else
		fftwf_execute(this->fwd);
#endif

		const MYReal norm_factor = 1./(gsize[0]*gsize[1]*gsize[2]);
		const MYReal h = grid.getBlocksInfo().front().h_gridpoint;

        _solveReproject((mycomplex *)this->data, gsize[0], gsize[1], gsize[2], gsize[2]/2+1, norm_factor, h);

#ifndef _FFTW_SP_COMP_
		fftw_execute(this->bwd);
#else
		fftwf_execute(this->bwd);
#endif

		this->_fftw2cub(this->data, grid, gsize[0], gsize[1], gsize[2], gsize[2]/2+1);

		//profiler.printSummary();
	}



};
