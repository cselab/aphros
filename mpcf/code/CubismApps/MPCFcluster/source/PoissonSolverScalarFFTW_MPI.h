#pragma once

#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include <algorithm>

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

using namespace std;

#include <BlockInfo.h>
#include <Profiler.h>

#include <Types.h>

// switch on initial energy spectrum taken from Comte-Bellot-Corrsin experiment for simulation
/* #define CBC_EXPERIMENT */

class FFTWBaseMPI
{
    static int registered_objects;
    static bool initialized; //a la singleton

    static void _setup(const int desired_threads)
    {
        if (!initialized)
        {
            initialized = true;

            int supported_threads;
            MPI_Query_thread(&supported_threads);

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

        registered_objects++;
    }

public:

    FFTWBaseMPI(const int desired_threads) { _setup(desired_threads); }

    static void dispose()
    {
        registered_objects--;

        if (registered_objects == 0)
        {
#ifndef _FFTW_SP_COMP_
            fftw_mpi_cleanup();
#else
            fftwf_mpi_cleanup();
#endif
        }
    }
};

template<typename TGrid, typename TStreamer>
class PoissonSolverScalarFFTW_MPI : FFTWBaseMPI
{
    Profiler profiler;

protected:
    typedef typename TGrid::BlockType BlockType;

    bool initialized;

    myplan fwd, bwd;
    ptrdiff_t alloc_local, local_n0, local_0_start, local_n1, local_1_start;
    MYReal * data; // rhs in _setup, out in cub2fftw and fftw2cub

protected:

    virtual void _setup(MYReal *& rhs , const size_t nx, const size_t ny, const size_t nz, MPI_Comm comm)
    {
        if (!initialized)
        {
            initialized = true;
#ifndef _FFTW_SP_COMP_
            alloc_local = fftw_mpi_local_size_3d_transposed(nx, ny, nz/2+1, comm, &local_n0, &local_0_start, &local_n1, &local_1_start);

            rhs = fftw_alloc_real(2*alloc_local);

            fwd = fftw_mpi_plan_dft_r2c_3d(nx, ny, nz, rhs, (mycomplex *)rhs, comm, FFTW_MPI_TRANSPOSED_OUT | FFTW_MEASURE);
            bwd = fftw_mpi_plan_dft_c2r_3d(nx, ny, nz, (mycomplex *)rhs, rhs, comm, FFTW_MPI_TRANSPOSED_IN | FFTW_MEASURE);
#else
            alloc_local = fftwf_mpi_local_size_3d_transposed(nx, ny, nz/2+1, comm, &local_n0, &local_0_start, &local_n1, &local_1_start);

            rhs = fftwf_alloc_real(2*alloc_local);

            fwd = fftwf_mpi_plan_dft_r2c_3d(nx, ny, nz, rhs, (mycomplex *)rhs, comm, FFTW_MPI_TRANSPOSED_OUT | FFTW_MEASURE);
            bwd = fftwf_mpi_plan_dft_c2r_3d(nx, ny, nz, (mycomplex *)rhs, rhs, comm, FFTW_MPI_TRANSPOSED_IN | FFTW_MEASURE);
#endif
        }
    }

    void _cub2fftw(TGrid& grid, MYReal * out, const size_t nx, const size_t ny, const size_t nz, const size_t nz_hat) const
    {
        vector<BlockInfo> local_infos = grid.getResidentBlocksInfo();

        const size_t N = local_infos.size();

        const size_t mybpd[3] = {static_cast<size_t>(grid.getResidentBlocksPerDimension(0)), static_cast<size_t>(grid.getResidentBlocksPerDimension(1)), static_cast<size_t>(grid.getResidentBlocksPerDimension(2))};

#pragma omp parallel for
        for(int i=0; i<N; ++i)
        {
            const BlockInfo info = local_infos[i];
            BlockType& b = *(BlockType*)local_infos[i].ptrBlock;

            const size_t offset = TStreamer::channels*(BlockType::sizeZ*info.index[2]+nz_hat*2*( BlockType::sizeY*info.index[1]+mybpd[1]*BlockType::sizeY*BlockType::sizeX*info.index[0]));


            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        const size_t dest_index = offset + TStreamer::channels*(iz + 2*nz_hat*( iy + BlockType::sizeY*mybpd[1]*ix ));

                        //assert(dest_index>=0 && dest_index<nx*ny*nz_hat*2); // linking error with openmp: http://lists.apple.com/archives/xcode-users/2011/Oct/msg00252.html
                        TStreamer::operate(b.data[iz][iy][ix], &out[dest_index]);
                    }
        }
    }


    /*-----------------------------------------------------------------------------------------------------
     * initialize velocity field for homogeneous isotropic turbulence in spectral space        Babak Hejazi
     *                                                                revised June 2015           rasthofer
     *-----------------------------------------------------------------------------------------------------*/
    void _solve_complex(mycomplex * in_out,       // Fourier coefficients of velocity field
            const size_t nx,          // number of points in x-direction (=number of grid cells in x-direction) (physical space)
            const size_t ny,          // number of points in y-direction (=number of grid cells in y-direction) (physical space)
            const size_t nz,          // number of points in z-direction (=number of grid cells in z-direction) (physical space)
            const size_t nz_hat,      // number of points in z-direction (spectral space): only one half of the spectral space need to be computed/stored,
            // the remaining part is obtained from the symmetry condition (u(k)=u*(k), since field is real in physical space)
            const MYReal norm_factor, // normalization of Fourier transformation by 1/(nz*ny*nz)
            const MYReal h            // cell size
            )
    {
        // subsequent implementations basically follow
        // R.S. Rogallo, Numerical experiments in homogeneous turbulence, NASA Tech. Memo. TM 81315, 1981.
        // S.S. Collis, Multiscale methods for turbulence simulation and control, Technical Report 034, MEMS, Rice University, 2002.
        // please check these references for further details

        srand48(M_PI);

        if (TStreamer::channels != 1)
        {
            cout << "PoissonSolverScalar::PoissonSolverScalar(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 1. Aborting\n";
            abort();
        }

        // usually the hit is solve in a 2pi-periodic domain
        // in this case, the wave factor would be zero
        // otherwise we scale to the appropriate length
        // TODO: (Ursula) validate
        const MYReal waveFactX = 2.0*M_PI/(nx*h);
        const MYReal waveFactY = 2.0*M_PI/(ny*h);
        const MYReal waveFactZ = 2.0*M_PI/(nz*h);

#pragma omp parallel for
        // loop all wave numbers
        for(int j = 0; j<local_n1; ++j) // local_n1=ny when using only one rank/mpi process
            for(int i=0; i<nx; ++i)
                for(int k = 0; k<nz_hat; ++k)
                {
                    // position in in_out vector
                    const size_t linidx = (j*nx+i)*nz_hat + k;

                    // transfer indices into wavenumber vector
                    const int k1 = (i <= nx/2) ? i : -(nx-i);
                    const int k2 = (local_1_start+j <= ny/2) ? local_1_start+j : -(ny-(local_1_start+j));
                    const int k3 = (k <= nz/2) ? k : -(nz-k);
                    // note once again that waveFacti=0 for a 2pi-periodic domain
                    const MYReal rk1 = k1*waveFactX;
                    const MYReal rk2 = k2*waveFactY;
                    const MYReal rk3 = k3*waveFactZ;

                    // get wavenumber
                    const MYReal K = sqrt(rk1*rk1+rk2*rk2+rk3*rk3);

                    // define energy spectrum
# ifndef CBC_EXPERIMENT
                    // wave number at which E(k) is maximal
                    /* const MYReal K0 = 4.0; */
                    const MYReal K0 = 10.0;
                    // initial rms velocity
                    /* const double u0_rms = 0.1; */
                    const double u0_rms = 0.0045;
                    // energy spectrum
                    const MYReal EK = K==0? 0 : 16.0 * sqrt(2.0/M_PI) * u0_rms * u0_rms/K0 * pow(K/K0,4.0) * exp(-2.0*pow(K/K0,2.0));
# else
                    // COMTE-BELLOT-CORRSIN EXPERIMENT

                    //----------------------------------------
                    // set-up wave numbers
                    //----------------------------------------

                    // wave number given from experiment [cm-1]
                    // have to be transferred to [m-1] and non-dimensionalized
                    std::vector<double> k_exp(19);
                    k_exp[0] = 0.20;
                    k_exp[1] = 0.25;
                    k_exp[2] = 0.30;
                    k_exp[3] = 0.40;
                    k_exp[4] = 0.50;
                    k_exp[5] = 0.70;
                    k_exp[6] = 1.00;
                    k_exp[7] = 1.50;
                    k_exp[8] = 2.00;
                    k_exp[9] = 2.50;
                    k_exp[10] = 3.00;
                    k_exp[11] = 4.00;
                    k_exp[12] = 6.00;
                    k_exp[13] = 8.00;
                    k_exp[14] = 10.00;
                    k_exp[15] = 12.50;
                    k_exp[16] = 15.00;
                    k_exp[17] = 17.50;
                    k_exp[18] = 20.00;

                    // non-dimensionalize wave number (according to Collis 2002)
                    // grid size of experiment
                    const double M = 0.0508;
                    // domain length
                    const double L = 10.0 * M;
                    // reference length
                    const double L_ref = L / (2.0 * M_PI);

                    for (std::size_t rr = 0; rr < k_exp.size(); rr ++)
                        k_exp[rr] *= (L_ref/0.01);

                    //----------------------------------------
                    // set-up energy
                    //----------------------------------------

                    // energy spectrum given from experiment [cm3/s2]
                    // have to be transferred to [m3/s2] and non-dimensionalized
                    std::vector<double> E_exp(19);
                    E_exp[0] = 129.0;
                    E_exp[1] = 230.0;
                    E_exp[2] = 322.0;
                    E_exp[3] = 435.0;
                    E_exp[4] = 457.0;
                    E_exp[5] = 380.0;
                    E_exp[6] = 270.0;
                    E_exp[7] = 168.0;
                    E_exp[8] = 120.0;
                    E_exp[9] = 89.0;
                    E_exp[10] = 70.3;
                    E_exp[11] = 47.0;
                    E_exp[12] = 24.7;
                    E_exp[13] = 12.6;
                    E_exp[14] = 7.42;
                    E_exp[15] = 3.96;
                    E_exp[16] = 2.33;
                    E_exp[17] = 1.34;
                    E_exp[18] = 0.80;

                    // non-dimensionalize energy spectrum
                    // inlet velocity of experiment
                    const double U_0 = 10.0;
                    // reference time
                    const double t_ref = 64.0 * M / U_0;

                    for (std::size_t rr = 0; rr < E_exp.size(); rr ++)
                        E_exp[rr] *= ((0.01*0.01*0.01)*(t_ref*t_ref)/(L_ref*L_ref*L_ref));

                    // the smallest value for k=1, the next smaller 1.41
                    // the following required extrapolation for k<1.6 yields
                    // negative values for k=1, which are not physical
                    // therefore, energy is set to zero for this case

                    // determine position of k
                    int position = -1;
                    for (std::size_t rr = 0; rr < k_exp.size(); rr ++)
                    {
                        if (K < k_exp[rr])
                        {
                            position = rr;
                            break;
                        }
                    }

                    // intialize energy
                    MYReal EK = 0.0;

                    if (position > 0)
                        // interpolate energy
                        EK = E_exp[position] + (E_exp[position-1] - E_exp[position]) / (k_exp[position] - k_exp[position-1]) * (k_exp[position] - K);
                    else if (position == 0)
                        // extrapolate energy
                        EK = E_exp[position+1] - (E_exp[position+1] - E_exp[position]) / (k_exp[position+1] - k_exp[position]) * (k_exp[position+1] - K);
                    else // position == -1: only possible for the present DNS which resolves the entire spectrum and even more
                        EK = 0.0; // TODO: (Ursula) is extrapolation better?

                    // see above (this corresponds to K=1 and K=0)
                    if (EK < 0.0)
                        EK = 0.0;

# endif

                    // get random phase angles between 0 and 2pi
                    const MYReal theta1 = drand48()*2.0*M_PI;
                    const MYReal theta2 = drand48()*2.0*M_PI;
                    const MYReal phi    = drand48()*2.0*M_PI;

                    const MYReal kfac = K==0? 0.0 : sqrt(EK/(2.0*M_PI*K*K));

                    const MYReal alpha_r = kfac*cos(theta1)*cos(phi);
                    const MYReal alpha_i = kfac*sin(theta1)*cos(phi);

                    const MYReal beta_r = kfac*cos(theta2)*sin(phi);
                    const MYReal beta_i = kfac*sin(theta2)*sin(phi);

                    const MYReal fac0 = K*sqrt(rk1*rk1+rk2*rk2);
                    const MYReal inv_fac0 = fac0==0 ? 0 : (MYReal)1/fac0;

                    // distinguish velocity components
                    const int comp = TStreamer::get_comp();

                    if (comp==0)
                    {
                        in_out[linidx][0] = K==0? 0 : inv_fac0 * (alpha_r * K * rk2 + beta_r * rk1 * rk3);
                        in_out[linidx][1] = K==0? 0 : inv_fac0 * (alpha_i * K * rk2 + beta_i * rk1 * rk3);
                    }
                    else if (comp==1)
                    {
                        in_out[linidx][0] = K==0? 0 : inv_fac0 * (-alpha_r * K * rk1 + beta_r * rk2 * rk3);
                        in_out[linidx][1] = K==0? 0 : inv_fac0 * (-alpha_i * K * rk1 + beta_i * rk2 * rk3);
                    }
                    else if (comp==2)
                    {
                        in_out[linidx][0] = K==0? 0 : -(beta_r*sqrt(k1*k1+k2*k2))/K;
                        in_out[linidx][1] = K==0? 0 : -(beta_i*sqrt(k1*k1+k2*k2))/K;
                    }
                    else
                    {
                        cout << "specify velocity component ..." << std::endl;
                        abort();
                    }
                }
    }

    virtual void _solve(mycomplex * in_out, const size_t nx, const size_t ny, const size_t nz, const size_t nz_hat, const MYReal norm_factor, const MYReal h)
    {
        if (TStreamer::channels != 1)
        {
            cout << "PoissonSolverScalar::PoissonSolverScalar(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 1. Aborting\n";
            abort();
        }

        const MYReal h2 = h*h;
        const MYReal factor = h2*norm_factor;

#pragma omp parallel for
        for(int j = 0; j<local_n1; ++j)
            for(int i=0; i<nx; ++i)
                for(int k = 0; k<nz_hat; ++k)
                {
                    const int linidx = (j*nx+i)*nz_hat + k;
                    //assert(linidx >=0 && linidx<nx*ny*nz_hat); // linking error with openmp

                    const MYReal denom = 32.*(cos(2.*M_PI*i/nx)+cos(2.*M_PI*(local_1_start+j)/ny)+
                            cos(2.0*M_PI*k/nz))-
                    2.*(cos(4.*M_PI*i/nx)+ cos(4.*M_PI*(local_1_start+j)/ny)+
                            cos(4.*M_PI*k/nz))-90.;

                    const MYReal inv_denom = (denom==0)? 0.:1./denom;
                    const MYReal fatfactor = 12. * inv_denom * factor;

                    in_out[linidx][0] *= fatfactor;
                    in_out[linidx][1] *= fatfactor;
                }

        //this is sparta!
        if (local_1_start == 0)
            in_out[0][0] = in_out[0][1] = 0;
    }

    virtual void _solveSpectral(mycomplex * in_out, const size_t nx, const size_t ny, const size_t nz, const size_t nz_hat, const MYReal norm_factor, const MYReal h)
    {
        if (TStreamer::channels != 1)
        {
            cout << "PoissonSolverScalar::PoissonSolverScalar(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 1. Aborting\n";
            abort();
        }

        const MYReal waveFactX = 2.0*M_PI/(nx*h);
        const MYReal waveFactY = 2.0*M_PI/(ny*h);
        const MYReal waveFactZ = 2.0*M_PI/(nz*h);

#pragma omp parallel for
        for(int j = 0; j<local_n1; ++j)
            for(int i=0; i<nx; ++i)
                for(int k = 0; k<nz_hat; ++k)
                {
                    const int linidx = (j*nx+i)*nz_hat + k;
                    //assert(linidx >=0 && linidx<nx*ny*nz_hat); // linking error with openmp

                    const int kx = (i <= nx/2) ? i : -(nx-i);
                    const int ky = (local_1_start+j <= ny/2) ? local_1_start+j : -(ny-(local_1_start+j));
                    const int kz = (k <= nz/2) ? k : -(nz-k);
                    const MYReal rkx = kx*waveFactX;
                    const MYReal rky = ky*waveFactY;
                    const MYReal rkz = kz*waveFactZ;

                    const MYReal kinv = (kx==0 && ky==0 && kz==0) ? 0.0 : -1.0/(rkx*rkx+rky*rky+rkz*rkz);
                    in_out[linidx][0] *= kinv*norm_factor;
                    in_out[linidx][1] *= kinv*norm_factor;
                }

        //this is sparta!
        if (local_1_start == 0)
            in_out[0][0] = in_out[0][1] = 0;
    }

    void _fftw2cub(MYReal * out, TGrid& grid, const size_t nx, const size_t ny, const size_t nz, const size_t nz_hat) const
    {
        vector<BlockInfo> local_infos = grid.getResidentBlocksInfo();

        const size_t N = local_infos.size();

        const size_t mybpd[3] = {static_cast<size_t>(grid.getResidentBlocksPerDimension(0)), static_cast<size_t>(grid.getResidentBlocksPerDimension(1)), static_cast<size_t>(grid.getResidentBlocksPerDimension(2))};

#pragma omp parallel for
        for(int i=0; i<N; ++i)
        {
            const BlockInfo info = local_infos[i];
            BlockType& b = *(BlockType*)local_infos[i].ptrBlock;

            const size_t offset = TStreamer::channels*(BlockType::sizeZ*info.index[2] + nz_hat*2*( BlockType::sizeY*info.index[1]+mybpd[1]*BlockType::sizeY*BlockType::sizeX*info.index[0]));

            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        const size_t src_index = offset + TStreamer::channels*(iz + 2*nz_hat*( iy + BlockType::sizeY*mybpd[1]*ix ));

                        //assert(src_index>=0 && src_index<nx*ny*nz_hat*2); // linking error with openmp
                        TStreamer::operate(&out[src_index], b.data[iz][iy][ix]);
                    }
        }
    }

    /*----------------------------------------------------------------------------------
     * compute energy spectrum from velocity field in spectral space        Babak Hejazi
     *                                             revised June 2015           rasthofer
     *----------------------------------------------------------------------------------*/
    void _fftw2cub_halfway(mycomplex * in_out,       // Fourier coefficients of velocity field
            const MPI_Comm comm,
            const size_t nx,          // number of points in x-direction (=number of grid cells in x-direction) (physical space)
            const size_t ny,          // number of points in y-direction (=number of grid cells in y-direction) (physical space)
            const size_t nz,          // number of points in z-direction (=number of grid cells in z-direction) (physical space)
            const size_t nz_hat,      // number of points in z-direction (spectral space): only one half of the spectral space need to be computed/stored,
            // the remaining part is obtained from the symmetry condition (u(k)=u*(k), since field is real in physical space)
            const MYReal norm_factor, // normalization of Fourier transformation by 1/(nz*ny*nz)
            const MYReal h,           // cell length
            const int step=0,
            const int comp=0)
    {
        // only one channel can be evaluated
        if (TStreamer::channels != 1)
        {
            cout << "PoissonSolverScalar::PoissonSolverScalar(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 1. Aborting\n";
            abort();
        }

        // usually the hit is solve in a 2pi-periodic domain
        // in this case, the wave factor would be zero
        // otherwise we scale to the appropriate length
        // TODO: (Ursula) validate and see also below at Nbins
        const MYReal waveFactX = 2.0*M_PI/(nx*h);
        const MYReal waveFactY = 2.0*M_PI/(ny*h);
        const MYReal waveFactZ = 2.0*M_PI/(nz*h);

        // note the maximal wavenumber occurring is sqrt(3)*nx/2 (assuming nx=ny=nz)
        // since bins of size one, ranging from k-0.5 to k+0.5, are used, we round up to the next integer
        if (nx!=ny or nx!=nz)
        {
            cout << "cubic domain assumed for hit" << endl;
            abort();
        }
        const MYReal binSize = 1.0;
        const int Nbins = ceil(sqrt(3)*nx/2.0*waveFactX)+1; // nx/2 rounded up plus one bin for k=0 (included waveFactX for scaling)
        MYReal energy[Nbins];
        for(int i=0; i<Nbins; i++)
            energy[i]=0.0;

#pragma omp parallel for
        // loop all wave numbers
        for(int j = 0; j<local_n1; ++j) // local_n1=ny when using only one rank/mpi process
            for(int i=0; i<nx; ++i)
                for(int k = 0; k<nz_hat; ++k)
                {
                    // position in in_out vector
                    const int linidx = (j*nx+i)*nz_hat + k;

                    // transfer indices into wavenumber vector
                    const int k1 = (i <= nx/2) ? i : -(nx-i);
                    const int k2 = (local_1_start+j <= ny/2) ? local_1_start+j : -(ny-(local_1_start+j)); // local_1_start=0 when using only one rank/mpi process
                    const int k3 = (k <= nz/2) ? k : -(nz-k);
                    // note once again that waveFacti=0 for a 2pi-periodic domain
                    const MYReal rk1 = k1*waveFactX;
                    const MYReal rk2 = k2*waveFactY;
                    const MYReal rk3 = k3*waveFactZ;

                    // get wavenumber
                    const MYReal K = sqrt(rk1*rk1+rk2*rk2+rk3*rk3);

                    // find corresponding bin
                    int binID = -1;
                    for (size_t r=0; r<Nbins; r++)
                    {
                        if (K<((MYReal) r + 0.5))
                        {
                            binID = r;
                            break;
                        }
                    }
                    if (binID == -1)
                    {
                        cout << "Could not find bin for wavenumber contribution!" << endl;
                        abort();
                    }

#pragma omp critical
                    // no factor 0.5 for energy since we perform evaluation only for one half of the Fourier space
                    // the second half lead to the same EK values since it is given by the conjugate
                    // complex values of the evaluated one (meaning EK * 2, i.e. 0.5*2=1)
                    energy[binID] += (pow(in_out[linidx][0],2)+pow(in_out[linidx][1],2))*norm_factor*norm_factor;
                }

        // sum up contributions of all ranks
        MYReal g_energy[Nbins];

        MPI_Reduce(&energy[0], &g_energy[0], Nbins, MPI_DOUBLE, MPI_SUM, 0, comm);

        // TODO: (Ursula) improve output
        // write energy spectrum to file
        // further compute the turbulent kinetic energy
        MYReal tke = 0.0;
        // the Taylor-micro scale
        MYReal lambda = 0.0;
        MYReal eps = 0.0;
        // as well as the Taylor-micro scale Reynolds number
        MYReal Re_lambda = 0.0;

        int myrank;
        MPI_Comm_rank(comm, &myrank);
        if (0 == myrank)
        {
            FILE * f_ke;

            char buf[256];
            sprintf(buf, "spectrum_%d_%d.dat", step, comp);
            f_ke = fopen(buf, "w");

            for(int i=0; i<Nbins; i++)
            {
                fprintf(f_ke, "%e %e\n", i*binSize, g_energy[i]);

#if 0
                // add to tke and dissipation eps
                tke += g_energy[i];
                eps += ((MYReal) i*binSize) * ((MYReal) i*binSize) g_energy[i];
#endif
            }

            fclose(f_ke);

            // TODO: (Ursula) this does not work this way since I need the contributions of the other velocity components
            //                shifted evaluation to python postprocessing script for the moment
# if 0
            FILE * f_int;

            char buf[256];
            sprintf(buf, "integral_spectrum_%d.dat", step, comp);
            f_ke = fopen(buf, "w");

            // get Taylor-micro scale
            lambda = sqrt(5.0*tke/eps);
            MYReal uprime = sqrt(2.0/3.0*tke);
            // Taylor-micro scale Reynolds number
            // TODO: (Ursula) get average density and viscosity
            Re_lambda = 1000.0 * uprime * lmbda / 0.001002;

            // write to file
            fprintf(f_int, "tke %e\n", tke);
#endif
        }
    }

public:

    PoissonSolverScalarFFTW_MPI(const int desired_threads): FFTWBaseMPI(desired_threads), initialized(false) {	}

    void solve(TGrid& grid, const bool spectral = false)
    {
        const int bs[3] = {BlockType::sizeX, BlockType::sizeY, BlockType::sizeZ};
        const size_t gsize[3] = {grid.getBlocksPerDimension(0)*bs[0], grid.getBlocksPerDimension(1)*bs[1], grid.getBlocksPerDimension(2)*bs[2]};

        MPI_Comm comm = grid.getCartComm();

        //profiler.push_start("SETUP");
        _setup(data, gsize[0], gsize[1], gsize[2], comm);
        //profiler.pop_stop();

        //profiler.push_start("CUB2FFTW");
        _cub2fftw(grid, data, local_n0, gsize[1], gsize[2], gsize[2]/2+1);
        //profiler.pop_stop();

        //profiler.push_start("FFTW FORWARD");
#ifndef _FFTW_SP_COMP_
        fftw_execute(fwd);
#else
        fftwf_execute(fwd);
#endif
        //profiler.pop_stop();

        const MYReal norm_factor = 1./(gsize[0]*gsize[1]*gsize[2]);
        const MYReal h = grid.getBlocksInfo().front().h_gridpoint;

        //profiler.push_start("SOLVE");
        if(spectral)
            _solveSpectral((mycomplex *)data, gsize[0], gsize[1], gsize[2], gsize[2]/2+1, norm_factor, h);
        else
            _solve((mycomplex *)data, gsize[0], gsize[1], gsize[2], gsize[2]/2+1, norm_factor, h);
        //profiler.pop_stop();

        //profiler.push_start("FFTW INVERSE");
#ifndef _FFTW_SP_COMP_
        fftw_execute(bwd);
#else
        fftwf_execute(bwd);
#endif
        //profiler.pop_stop();

        //profiler.push_start("FFTW2CUB");
        _fftw2cub(data, grid, gsize[0], gsize[1], gsize[2], gsize[2]/2+1);
        //profiler.pop_stop();

        //profiler.printSummary();
    }

    void solve_halfway_back(TGrid& grid, const bool spectral = false)
    {
        const int bs[3] = {BlockType::sizeX, BlockType::sizeY, BlockType::sizeZ};
        const size_t gsize[3] = {
            static_cast<size_t>(grid.getBlocksPerDimension(0)*bs[0]),
            static_cast<size_t>(grid.getBlocksPerDimension(1)*bs[1]),
            static_cast<size_t>(grid.getBlocksPerDimension(2)*bs[2])};

        MPI_Comm comm = grid.getCartComm();

        profiler.push_start("SETUP");
        _setup(data, gsize[0], gsize[1], gsize[2], comm);
        profiler.pop_stop();

        profiler.push_start("CUB2FFTW");
        _cub2fftw(grid, data, local_n0, gsize[1], gsize[2], gsize[2]/2+1);
        profiler.pop_stop();

        const MYReal norm_factor = 1./(gsize[0]*gsize[1]*gsize[2]);
        const MYReal h = grid.getBlocksInfo().front().h_gridpoint;

        profiler.push_start("SOLVE");
        _solve_complex((mycomplex *)data, gsize[0], gsize[1], gsize[2], gsize[2]/2+1, norm_factor, h);
        profiler.pop_stop();

        profiler.push_start("FFTW INVERSE");
#ifndef _FFTW_SP_COMP_
        fftw_execute(bwd);
#else
        fftwf_execute(bwd);
#endif
        profiler.pop_stop();

        profiler.push_start("FFTW2CUB");
        _fftw2cub(data, grid, gsize[0], gsize[1], gsize[2], gsize[2]/2+1);
        profiler.pop_stop();

        //profiler.printSummary();
    }

    void solve_halfway_forward(TGrid& grid, const bool bPrint = false)
    {
        const int bs[3] = {BlockType::sizeX, BlockType::sizeY, BlockType::sizeZ};
        const size_t gsize[3] = {static_cast<size_t>(grid.getBlocksPerDimension(0)*bs[0]), static_cast<size_t>(grid.getBlocksPerDimension(1)*bs[1]), static_cast<size_t>(grid.getBlocksPerDimension(2)*bs[2])};

        MPI_Comm comm = grid.getCartComm();

        profiler.push_start("SETUP");
        _setup(data, gsize[0], gsize[1], gsize[2], comm);
        profiler.pop_stop();

        profiler.push_start("CUB2FFTW");
        _cub2fftw(grid, data, local_n0, gsize[1], gsize[2], gsize[2]/2+1);
        profiler.pop_stop();

        profiler.push_start("FFTW INVERSE");
#ifndef _FFTW_SP_COMP_
        fftw_execute(fwd);
#else
        fftwf_execute(fwd);
#endif
        profiler.pop_stop();

        const MYReal norm_factor = 1./(gsize[0]*gsize[1]*gsize[2]);
        const MYReal h = grid.getBlocksInfo().front().h_gridpoint;

        profiler.push_start("FFTW2CUB HALFWAY");
        if (bPrint) _fftw2cub_halfway((mycomplex *)data, comm, gsize[0], gsize[1], gsize[2], gsize[2]/2+1, norm_factor, h);
        profiler.pop_stop();

        profiler.push_start("FFTW2CUB");
        _fftw2cub(data, grid, gsize[0], gsize[1], gsize[2], gsize[2]/2+1);
        profiler.pop_stop();

        //profiler.printSummary();
    }

    void calc_spectra(TGrid& grid, const int step, const int comp)
    {
        const int bs[3] = {BlockType::sizeX, BlockType::sizeY, BlockType::sizeZ};
        const size_t gsize[3] = {static_cast<size_t>(grid.getBlocksPerDimension(0)*bs[0]), static_cast<size_t>(grid.getBlocksPerDimension(1)*bs[1]), static_cast<size_t>(grid.getBlocksPerDimension(2)*bs[2])};

        MPI_Comm comm = grid.getCartComm();

        profiler.push_start("SETUP");
        _setup(data, gsize[0], gsize[1], gsize[2], comm);
        profiler.pop_stop();

        profiler.push_start("CUB2FFTW");
        // TODO: (Ursula) Why here local_n0 and in _fftw2cub_halfway local_n1?
        _cub2fftw(grid, data, local_n0, gsize[1], gsize[2], gsize[2]/2+1);
        profiler.pop_stop();

        profiler.push_start("FFTW INVERSE");
#ifndef _FFTW_SP_COMP_
        fftw_execute(fwd);
#else
        fftwf_execute(fwd);
#endif
        profiler.pop_stop();

        const MYReal norm_factor = 1./(gsize[0]*gsize[1]*gsize[2]);
        const MYReal h = grid.getBlocksInfo().front().h_gridpoint;

        profiler.push_start("FFTW2CUB HALFWAY");
        _fftw2cub_halfway((mycomplex *)data, comm, gsize[0], gsize[1], gsize[2], gsize[2]/2+1, norm_factor, h, step, comp);
        profiler.pop_stop();

        //profiler.printSummary();
    }

    void dispose()
    {
        if (initialized)
        {
            initialized = false;

#ifndef _FFTW_SP_COMP_
            fftw_destroy_plan(fwd);
            fftw_destroy_plan(bwd);

            fftw_free(data);
#else
            fftwf_destroy_plan(fwd);
            fftwf_destroy_plan(bwd);

            fftwf_free(data);
#endif
            FFTWBaseMPI::dispose();
        }
    }
};
