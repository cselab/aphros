/* File:   evalHit.cpp */
/* Date:   July 2015*/
/* Author: Ursula Rasthofer */
/* Tag:    evaluatiosn for hit */
/* Copyright 2015 ETH Zurich. All Rights Reserved. */

#include "evalHit.h"
#include "../HDF5Reader_MPI.h"

#include <omp.h>

/*---------------------------------------------------------------------------*
 | constructor                                           rasthofer July 2015
 *---------------------------------------------------------------------------*/
evalHit::evalHit(ArgumentParser& myparser):
evalBase(),
FFTWVectorMPI(omp_get_max_threads()),
myparser_(myparser)
{
  return;
}


/*---------------------------------------------------------------------------*
 | evaluate solution in spectral space                   rasthofer July 2015
 *---------------------------------------------------------------------------*/
void evalHit::compute()
{
  //---------------------------------------
  // general setup
  //---------------------------------------

  // get grid size
  myparser_.set_strict_mode();
  const std::size_t nx = myparser_("NX").asInt();
  const std::size_t ny = myparser_("NY").asInt();
  const std::size_t nz = myparser_("NZ").asInt();
  const MYReal maxextent = (MYReal)myparser_("extent").asDouble();
  myparser_.unset_strict_mode();

  const std::size_t nz_hat = nz/2+1;

  // cell size
  const MYReal h = maxextent / std::max( nx, std::max( ny, nz ) );
  // and normalization
  const MYReal norm_factor = 1.0/(nx*ny*nz);

  // get hdf channels for evaluations
  std::string hdfchannels = myparser_("-hdfchannels").asString("1");

  //--------------------------------------
  // setup fft solver
  //--------------------------------------

  // we always want to transfer the velocity
  _setup(data, nx, ny, nz);

  //--------------------------------------
  // get data
  //--------------------------------------

  std::size_t nchannels = 3;
  myparser_.set_strict_mode();
  std::string filename = myparser_("-vel").asString();
  myparser_.unset_strict_mode();

  if (hdfchannels.find('a') != std::string::npos)
    nchannels = 7;

  data_2 = new Real[local_n0 * ny * nz * nchannels];

  ReadHDF5_MPI(data_2, nx, ny, nz, local_n0, ny, nz, nchannels, filename);

  //---------------------------------------
  // do velocity
  //---------------------------------------

  // fill data
  if (nchannels == 3)
  {
    // transfer input to fftw structure
    // with padding (see documentation)
    for (int k = 0; k<nz; k++)
      for (int j = 0; j<ny; j++)
        for (int i = 0; i<local_n0; i++)
        {
          int index = nchannels * (i + local_n0 * (j + ny * k));
          int index_d = nchannels * (k + 2 * nz_hat * (j + ny *i));
          for (int idim = 0; idim<nchannels; idim++)
            data[index_d+idim] = (MYReal)(data_2[index+idim]);
        }

  }
  else
  {
    // TODO:
    abort();
  }

  // perform fftw
  fftw_execute(fwd);

  // compute spectrum

  // usually the hit is solve in a 2pi-periodic domain
  // in this case, the wave factor would be 1
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
  const int Nbins = ceil(sqrt(3.0)*nx/2.0*waveFactX)+1.0; // nx/2 rounded up plus one bin for k=0 (included waveFactX for scaling)
  MYReal energy[Nbins];
  MYReal energy_x[Nbins];
  MYReal energy_y[Nbins];
  MYReal energy_z[Nbins];
  MYReal energy_dil[Nbins];
  MYReal energy_sol[Nbins];
  MYReal energy_pp[Nbins];
  MYReal energy_rr[Nbins];
  for(int i=0; i<Nbins; i++)
  {
    energy[i]=0.0;
    energy_x[i]=0.0;
    energy_y[i]=0.0;
    energy_z[i]=0.0;
    energy_dil[i]=0.0;
    energy_sol[i]=0.0;
    energy_pp[i]=0.0;
    energy_rr[i]=0.0;
  }

#pragma omp parallel for
  // loop all wave numbers
  for(int j = 0; j<local_n1; ++j) // local_n1=ny when using only one rank/mpi process
    for(int i=0; i<nx; ++i)
      for(int k = 0; k<nz_hat; ++k)
      {
        // position in in_out vector
        const int linidx = 3 * ((j*nx+i)*nz_hat + k);

        // transfer indices into wavenumber vector
        const int k1 = (i <= nx/2) ? i : -(nx-i);
        const int k2 = (local_1_start+j <= ny/2) ? local_1_start+j : -(ny-(local_1_start+j)); // local_1_start=0 when using only one rank/mpi process
        const int k3 = (k <= nz/2) ? k : -(nz-k);
        // note once again that waveFacti=1 for a 2pi-periodic domain
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

        // no factor 0.5 for energy since we perform evaluation only for one half of the Fourier space
        // the second half lead to the same EK values since it is given by the conjugate
        // complex values of the evaluated one (meaning EK * 2, i.e. 0.5*2=1)
        const MYReal squx = (pow(((mycomplex *)data)[linidx][0],2.0)+pow(((mycomplex *)data)[linidx][1],2.0));
        const MYReal squy = (pow(((mycomplex *)data)[linidx+1][0],2.0)+pow(((mycomplex *)data)[linidx+1][1],2.0));
        const MYReal squz = (pow(((mycomplex *)data)[linidx+2][0],2.0)+pow(((mycomplex *)data)[linidx+2][1],2.0));

        // compute solenoidal and dilatational velocity parts and their spectra
        const MYReal divu_re = (rk1*((mycomplex *)data)[linidx][0] + rk2*((mycomplex *)data)[linidx+1][0] + rk3*((mycomplex *)data)[linidx+2][0]) * norm_factor;
        const MYReal divu_im = (rk1*((mycomplex *)data)[linidx][1] + rk2*((mycomplex *)data)[linidx+1][1] + rk3*((mycomplex *)data)[linidx+2][1]) * norm_factor;

        const MYReal ux_dil_re = (K>0.0) ? rk1/(K*K) * (divu_re * ((mycomplex *)data)[linidx][0] - divu_im * ((mycomplex *)data)[linidx][1]) * norm_factor : 0.0;
        const MYReal ux_dil_im = (K>0.0) ? rk1/(K*K) * (divu_re * ((mycomplex *)data)[linidx][1] + divu_im * ((mycomplex *)data)[linidx][0]) * norm_factor : 0.0;
        const MYReal uy_dil_re = (K>0.0) ? rk2/(K*K) * (divu_re * ((mycomplex *)data)[linidx+1][0] - divu_im * ((mycomplex *)data)[linidx+1][1]) * norm_factor : 0.0;
        const MYReal uy_dil_im = (K>0.0) ? rk2/(K*K) * (divu_re * ((mycomplex *)data)[linidx+1][1] + divu_im * ((mycomplex *)data)[linidx+1][0]) * norm_factor : 0.0;
        const MYReal uz_dil_re = (K>0.0) ? rk3/(K*K) * (divu_re * ((mycomplex *)data)[linidx+2][0] - divu_im * ((mycomplex *)data)[linidx+2][1]) * norm_factor : 0.0;
        const MYReal uz_dil_im = (K>0.0) ? rk3/(K*K) * (divu_re * ((mycomplex *)data)[linidx+2][1] + divu_im * ((mycomplex *)data)[linidx+2][0]) * norm_factor : 0.0;

        const MYReal ux_sol_re = ((mycomplex *)data)[linidx][0] * norm_factor - ux_dil_re;
        const MYReal ux_sol_im = ((mycomplex *)data)[linidx][1] * norm_factor - ux_dil_im;
        const MYReal uy_sol_re = ((mycomplex *)data)[linidx+1][0] * norm_factor - uy_dil_re;
        const MYReal uy_sol_im = ((mycomplex *)data)[linidx+1][1] * norm_factor - uy_dil_im;
        const MYReal uz_sol_re = ((mycomplex *)data)[linidx+2][0] * norm_factor - uz_dil_re;
        const MYReal uz_sol_im = ((mycomplex *)data)[linidx+2][1] * norm_factor - uz_dil_im;
#pragma omp critical
{
        energy[binID] += (squx+squy+squz)*norm_factor*norm_factor;
        energy_x[binID] += squx*norm_factor*norm_factor;
        energy_y[binID] += squy*norm_factor*norm_factor;
        energy_z[binID] += squz*norm_factor*norm_factor;
        energy_dil[binID] += (ux_dil_re*ux_dil_re + ux_dil_im*ux_dil_im + uy_dil_re*uy_dil_re + uy_dil_im*uy_dil_im + uz_dil_re*uz_dil_re + uz_dil_im*uz_dil_im);
        energy_sol[binID] += (ux_sol_re*ux_sol_re + ux_sol_im*ux_sol_im + uy_sol_re*uy_sol_re + uy_sol_im*uy_sol_im + uz_sol_re*uz_sol_re + uz_sol_im*uz_sol_im);
}

      }

  // clean up
  dispose();

  // do pressure
  if (hdfchannels.find('4') != std::string::npos)
    nchannels = 1;

  if (nchannels == 1 or nchannels == 7)
  {
    //--------------------------------
    // set up scalar fftw
    //--------------------------------

    FFTWScalarMPI::_setup(data, nx, ny, nz);

    //--------------------------------
    // get data
    //--------------------------------

    if (nchannels == 1)
    {
      myparser_.set_strict_mode();
      std::string filename = myparser_("-pres").asString();
      myparser_.unset_strict_mode();
    }

    // fill data if not all field have been read from one file
    if (nchannels == 1)
    {
      data_2 = new Real[local_n0 * ny * nz * nchannels];

      ReadHDF5_MPI(data_2, nx, ny, nz, local_n0, ny, nz, nchannels, filename);

      // transfer input to fftw structure
      // with padding (see documentation)
      for (int k = 0; k<nz; k++)
        for (int j = 0; j<ny; j++)
          for (int i = 0; i<local_n0; i++)
          {
            int index = nchannels * (i + local_n0 * (j + ny * k));
            int index_d = nchannels * (k + 2 * nz_hat * (j + ny *i));
            data[index_d] = (MYReal)(data_2[index]);
          }

    }
    else
    {
      // TODO:
      abort();
    }

    // perform fftw
    fftw_execute(fwd);

#pragma omp parallel for
    // loop all wave numbers
    for(int j = 0; j<local_n1; ++j) // local_n1=ny when using only one rank/mpi process
      for(int i=0; i<nx; ++i)
        for(int k = 0; k<nz_hat; ++k)
        {
          // position in in_out vector
          const int linidx = ((j*nx+i)*nz_hat + k);

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

          // get square of pressure
          const MYReal sqpp = (pow(((mycomplex *)data)[linidx][0],2.0)+pow(((mycomplex *)data)[linidx][1],2.0));
          // compute spectrum
#pragma omp critical
{
          energy_pp[binID] += sqpp*norm_factor*norm_factor;
}

        }
     // reset
     if (nchannels == 1) nchannels = 3;
     dispose();
  }

  // do density
  if (hdfchannels.find('0') != std::string::npos)
    nchannels = 1;

  if (nchannels == 1 or nchannels == 7)
  {
    //--------------------------------
    // set up scalar fftw
    //--------------------------------

    FFTWScalarMPI::_setup(data, nx, ny, nz);

    //--------------------------------
    // get data
    //--------------------------------

    if (nchannels == 1)
    {
      myparser_.set_strict_mode();
      std::string filename = myparser_("-density").asString();
      myparser_.unset_strict_mode();
    }

    // fill data if not all field have been read from one file
    if (nchannels == 1)
    {
      data_2 = new Real[local_n0 * ny * nz * nchannels];

      ReadHDF5_MPI(data_2, nx, ny, nz, local_n0, ny, nz, nchannels, filename);

      // transfer input to fftw structure
      // with padding (see documentation)
      for (int k = 0; k<nz; k++)
        for (int j = 0; j<ny; j++)
          for (int i = 0; i<local_n0; i++)
          {
            int index = nchannels * (i + local_n0 * (j + ny * k));
            int index_d = nchannels * (k + 2 * nz_hat * (j + ny *i));
            data[index_d] = (MYReal)(data_2[index]);
          }
    }
    else
    {
      // TODO:
      abort();
    }

    // perform fftw
    fftw_execute(fwd);

#pragma omp parallel for
    // loop all wave numbers
    for(int j = 0; j<local_n1; ++j) // local_n1=ny when using only one rank/mpi process
      for(int i=0; i<nx; ++i)
        for(int k = 0; k<nz_hat; ++k)
        {
          // position in in_out vector
          const int linidx = ((j*nx+i)*nz_hat + k);

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

          // get square of density
          const MYReal sqrr = (pow(((mycomplex *)data)[linidx][0],2.0)+pow(((mycomplex *)data)[linidx][1],2.0));
          // compute spectrum
#pragma omp critical
{
          energy_rr[binID] += sqrr*norm_factor*norm_factor;
}

        }
     // reset
     dispose();
  }

  //---------------------------------------
  // write to file
  //---------------------------------------

  // sum up contributions of all ranks
  MYReal g_energy[Nbins];
  MYReal g_energy_x[Nbins];
  MYReal g_energy_y[Nbins];
  MYReal g_energy_z[Nbins];
  MYReal g_energy_dil[Nbins];
  MYReal g_energy_sol[Nbins];
  MYReal g_energy_pp[Nbins];
  MYReal g_energy_rr[Nbins];

#ifndef _FFTW_SP_COMP_
  MPI::COMM_WORLD.Reduce(&energy[0], &g_energy[0], Nbins, MPI_DOUBLE, MPI::SUM, 0);
  MPI::COMM_WORLD.Reduce(&energy_x[0], &g_energy_x[0], Nbins, MPI_DOUBLE, MPI::SUM, 0);
  MPI::COMM_WORLD.Reduce(&energy_y[0], &g_energy_y[0], Nbins, MPI_DOUBLE, MPI::SUM, 0);
  MPI::COMM_WORLD.Reduce(&energy_z[0], &g_energy_z[0], Nbins, MPI_DOUBLE, MPI::SUM, 0);
  MPI::COMM_WORLD.Reduce(&energy_dil[0], &g_energy_dil[0], Nbins, MPI_DOUBLE, MPI::SUM, 0);
  MPI::COMM_WORLD.Reduce(&energy_sol[0], &g_energy_sol[0], Nbins, MPI_DOUBLE, MPI::SUM, 0);
  MPI::COMM_WORLD.Reduce(&energy_pp[0], &g_energy_pp[0], Nbins, MPI_DOUBLE, MPI::SUM, 0);
  MPI::COMM_WORLD.Reduce(&energy_rr[0], &g_energy_rr[0], Nbins, MPI_DOUBLE, MPI::SUM, 0);
#else
  MPI::COMM_WORLD.Reduce(&energy[0], &g_energy[0], Nbins, MPI_FLOAT, MPI::SUM, 0);
  MPI::COMM_WORLD.Reduce(&energy_x[0], &g_energy_x[0], Nbins, MPI_FLOAT, MPI::SUM, 0);
  MPI::COMM_WORLD.Reduce(&energy_y[0], &g_energy_y[0], Nbins, MPI_FLOAT, MPI::SUM, 0);
  MPI::COMM_WORLD.Reduce(&energy_z[0], &g_energy_z[0], Nbins, MPI_FLOAT, MPI::SUM, 0);
  MPI::COMM_WORLD.Reduce(&energy_dil[0], &g_energy_dil[0], Nbins, MPI_FLOAT, MPI::SUM, 0);
  MPI::COMM_WORLD.Reduce(&energy_sol[0], &g_energy_sol[0], Nbins, MPI_FLOAT, MPI::SUM, 0);
  MPI::COMM_WORLD.Reduce(&energy_pp[0], &g_energy_pp[0], Nbins, MPI_FLOAT, MPI::SUM, 0);
  MPI::COMM_WORLD.Reduce(&energy_rr[0], &g_energy_rr[0], Nbins, MPI_FLOAT, MPI::SUM, 0);
#endif

  if (MPI::COMM_WORLD.Get_rank()==0)
  {
    // write spectra to file
    FILE * f_ke;

    // get step for file name
    myparser_.set_strict_mode();
    std::string stepfile = myparser_("-step").asString();
    myparser_.unset_strict_mode();

    // further compute the turbulent kinetic energy
    MYReal tke = 0.0;
    // the Taylor-micro scale
    MYReal lambda = 0.0;
    MYReal eps = 0.0;
    // as well as the Taylor-micro scale Reynolds number
    MYReal Re_lambda = 0.0;

    // and the turbulent kinetic energies associated with the solenoidal and dilational velocity field
    MYReal Kss = 0.0;
    MYReal Kdd = 0.0;
    MYReal Chi = 0.0;

    char buf[256];
    std::string name_spec = "spectrum-"+stepfile+".dat";
    strcpy(&(buf[0]), name_spec.c_str());
    f_ke = fopen(buf, "w");

    fprintf(f_ke, "#    k            E_x            E_y            E_z            E            E_dil            E_sol            E_pp            E_rr\n");
    for(int i=0; i<Nbins; i++)
    {
      fprintf(f_ke, "%e   %e   %e   %e   %e   %e   %e   %e   %e\n", i*binSize, g_energy_x[i], g_energy_y[i], g_energy_z[i], g_energy[i], g_energy_dil[i], g_energy_sol[i], g_energy_pp[i], g_energy_rr[i]);

      // add to tke and dissipation eps
      tke += g_energy[i];
      eps += ((MYReal) i*binSize) * ((MYReal) i*binSize) * g_energy[i];
      Kss += g_energy_sol[i];
      Kdd += g_energy_dil[i];
    }

    fclose(f_ke);

    FILE * f_int;

    char buf2[256];
    name_spec = "integral_spectrum-"+stepfile+".dat";
    strcpy(&(buf2[0]), name_spec.c_str());
    f_int = fopen(buf2, "w");

    // get Taylor-micro scale
    lambda = sqrt(5.0*tke/eps);
    MYReal uprime = sqrt(2.0/3.0*tke);

    // get Chi
    Chi = Kdd / (Kdd + Kss);

    // Taylor-micro scale Reynolds number
    myparser_.set_strict_mode();
    const MYReal dens = (MYReal)myparser_("-rho1").asDouble();
    const MYReal visc = (MYReal)myparser_("-mu1").asDouble();
    myparser_.unset_strict_mode();
    Re_lambda = dens * uprime * lambda / visc;

    // write to file
    fprintf(f_int, "tke: %e\n", tke);
    fprintf(f_int, "eps: %e\n", eps);
    fprintf(f_int, "lambda: %e\n", lambda);
    fprintf(f_int, "uprime: %e\n", uprime);
    fprintf(f_int, "Re_lambda: %e\n", Re_lambda);
    fprintf(f_int, "Kss: %e\n", Kss);
    fprintf(f_int, "Kdd: %e\n", Kdd);
    fprintf(f_int, "Chi: %e\n", Chi);

    fclose(f_int);
  }

  return;
}


void evalHit::dispose()
{
  delete [] data_2;

  FFTWScalarMPI::dispose();

  return;
}
