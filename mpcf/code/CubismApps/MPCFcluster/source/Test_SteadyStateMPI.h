/*
 *  Test_SteadyStateMPI.h
 *  MPCFcluster
 *
 *  Created by Diego Rossinelli on 11/25/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#ifndef TEST_STEADYSTATEMPI_H_MT52JRKU
#define TEST_STEADYSTATEMPI_H_MT52JRKU


#include <cassert>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <ctime>

#include <mpi.h>

#include "Tests.h"
#include "BlockInfo.h"
#include "Types.h"
#include "GridMPI.h"
#include "BlockLabMPI.h"
#include "OutputProcessing.h"
#include "OutputProcessingMPI.h"


template <typename TGrid, typename TStepper, template <typename> class TSlice=SliceMPI>
class Test_SteadyStateMPI : public Simulation
{
  public:

    Test_SteadyStateMPI(const MPI_Comm comm, ArgumentParser& P) :
      grid(NULL), stepper(NULL), dumper(NULL),
      isroot(true),
      parser(P), restart_id(0), t(0.0), step_id(0),
      m_comm_world(comm)
  {
    int rank;
    MPI_Comm_rank(m_comm_world, &rank);
    this->isroot = (0 == rank);
    if (!this->isroot) this->VERBOSITY = 0;
  }
    virtual ~Test_SteadyStateMPI() {
      if (this->dumper)  delete dumper;
      if (this->stepper) delete stepper;
      if (this->grid)    delete grid;
    }

    inline void verbosity(const int v) { VERBOSITY = v; }

    virtual void setup();


  private:
    TGrid * grid;
    TStepper * stepper;
    OutputProcessing<TGrid,TSlice> * dumper;
    bool isroot;

    int BPDX, BPDY, BPDZ;
    int NSTEPS, SAVEPERIOD, VERBOSITY, REPORT_FREQ, REFRESHPERIOD;
    int LASTSAVE;
    Real TEND;

    int restart_id;
    bool BC_PERIODIC[3];
    Real t, dt;
    int step_id;

    ArgumentParser& parser;

    MPI_Comm m_comm_world;

    int XPESIZE, YPESIZE, ZPESIZE;

    virtual void _setEnvironment()
    {
      //parser.read_runtime_environment();

      VERBOSITY      = parser("verbosity").asInt();
      REPORT_FREQ    = parser("report").asInt();

      TEND           = parser("tend").asDouble();

      // OutputProcessing
      parser.set_strict_mode();
      dumper->m_dumpperiod   = parser("dumpperiod").asInt();
      dumper->m_dumpdt       = parser("dumpdt").asDouble();
      dumper->m_bIO          = parser("io").asBool();
      dumper->m_bVP          = parser("vp").asBool();
      dumper->m_bHDF         = parser("hdf").asBool();
      dumper->m_bHDF_SLICE   = parser("hdf_slice").asBool();
      dumper->m_heavySkipStep= parser("heavyskipstep").asInt();
      dumper->m_channels     = parser("channels").asString();
      parser.unset_strict_mode();
    }

    virtual void _setup_parameter();

    virtual void _ic();

    virtual void run();

    // case specific
    virtual void _print_case_header()
    {
      printf("////////////////////////////////////////////////////////////\n");
      printf("////////////      TEST STEADY STATE   MPI    ///////////////\n");
      printf("////////////////////////////////////////////////////////////\n");
      typedef typename TGrid::BlockType B;
      std::cout << "Domain size:   [" << this->BPDX*B::sizeX*XPESIZE;
      std::cout << " x " << this->BPDY*B::sizeY*YPESIZE;
      std::cout << " x " <<  this->BPDZ*B::sizeZ*ZPESIZE << "]" << std::endl;

      std::cout << "Domain extent: [" << Simulation_Environment::extents[0];
      std::cout << " x " << Simulation_Environment::extents[1];
      std::cout << " x " <<  Simulation_Environment::extents[2] << "]" << std::endl;
    }

    virtual void _init()
    {
      MPI_Barrier(m_comm_world);
    }
};


// class implementation
  template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SteadyStateMPI<TGrid,TStepper,TSlice>::_setup_parameter()
{
  parser.mute();

  parser.set_strict_mode();
  BPDX       = parser("-bpdx").asInt();
  BPDY           = parser("-bpdy").asInt(BPDX);
  BPDZ           = parser("-bpdz").asInt(BPDX);
  TEND       = parser("-tend").asDouble();
  parser.unset_strict_mode();

  // defaults

  VERBOSITY      = parser("-verbosity").asInt(0);
  Simulation_Environment::extent = parser("-extent").asDouble(1.0);

  // some post computations
  {
    Lab dummy;

    BC_PERIODIC[0] = dummy.is_xperiodic();
    BC_PERIODIC[1] = dummy.is_yperiodic();
    BC_PERIODIC[2] = dummy.is_zperiodic();

    Simulation_Environment::BC_PERIODIC[0] = BC_PERIODIC[0];
    Simulation_Environment::BC_PERIODIC[1] = BC_PERIODIC[1];
    Simulation_Environment::BC_PERIODIC[2] = BC_PERIODIC[2];
  }

  const int BPD_MAX = max(max(BPDX, BPDY), BPDZ);
  Simulation_Environment::extents[0] = Simulation_Environment::extent * BPDX / (double) BPD_MAX;
  Simulation_Environment::extents[1] = Simulation_Environment::extent * BPDY / (double) BPD_MAX;
  Simulation_Environment::extents[2] = Simulation_Environment::extent * BPDZ / (double) BPD_MAX;

  // some checks
  assert(TEND >= 0.0);
  assert(BPDX >= 1);
  assert(BPDY >= 1);
  assert(BPDZ >= 1);

  XPESIZE = this->parser("-xpesize").asInt(2);
  YPESIZE = this->parser("-ypesize").asInt(2);
  ZPESIZE = this->parser("-zpesize").asInt(2);

  const int bpdx = this->BPDX;
  const int bpdy = this->BPDY;
  const int bpdz = this->BPDZ;
  const int BPD_PE_MAX = 
    std::max(std::max(bpdx*XPESIZE, bpdy*YPESIZE), bpdz*ZPESIZE);
  Simulation_Environment::extents[0] = 
    Simulation_Environment::extent * 
    (bpdx*XPESIZE) / static_cast<double>(BPD_PE_MAX);
  Simulation_Environment::extents[1] = 
    Simulation_Environment::extent * 
    (bpdy*YPESIZE) / static_cast<double>(BPD_PE_MAX);
  Simulation_Environment::extents[2] = 
    Simulation_Environment::extent * 
    (bpdz*ZPESIZE) / static_cast<double>(BPD_PE_MAX);
}


  template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SteadyStateMPI<TGrid,TStepper,TSlice>::setup()
{
  _setup_parameter();

  if (this->isroot)
    _print_case_header();

  this->grid = new TGrid(
      XPESIZE, YPESIZE, ZPESIZE, 
      this->BPDX, this->BPDY, this->BPDZ, 
      Simulation_Environment::extent, m_comm_world);
  assert(this->grid != NULL);

  this->stepper = new TStepper(*(this->grid), this->parser, this->VERBOSITY);
  assert(this->stepper != NULL);

  this->dumper = new OutputProcessingMPI<TGrid,TSlice>(this->parser, *(this->grid), this->isroot);
  assert(this->dumper != NULL);
  this->dumper->register_all(*(this->grid));

  const std::string path = this->parser("-fpath").asString(".");
  this->_ic();
}

  template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SteadyStateMPI<TGrid,TStepper,TSlice>::_ic()
{
  std::vector<BlockInfo> vInfo = grid->getBlocksInfo();

  const double a2   = parser("a2").asDouble(0.0);

  typedef typename TGrid::BlockType TBlock;

#pragma omp parallel for
  for(int i=0; i<(int)vInfo.size(); i++)
  {
    BlockInfo info = vInfo[i];
    TBlock& b = *(TBlock*)info.ptrBlock;
    double* o = info.origin;

    for(int iz=0; iz<TBlock::sizeZ; iz++)
      for(int iy=0; iy<TBlock::sizeY; iy++)
        for(int ix=0; ix<TBlock::sizeX; ix++)
        {
          typedef const Real CReal;
          CReal e = Simulation_Environment::extent;
          Real x = o[0] + info.h * ix / TBlock::sizeX;
          Real y = o[1] + info.h * iy / TBlock::sizeY;
          x /= e;
          y /= e;
          b(ix, iy, iz).alpha2 = 0.5 * (1. + std::sin(x * 10 + y * y * 5));
        }
  }
}


  template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SteadyStateMPI<TGrid,TStepper,TSlice>::run()
{
  MPI_Barrier(m_comm_world);

  dt = parser("dt").asDouble(TEND / 100);
  int stepend = parser("stepend").asInt(100);

  while (true)
  {
    _setEnvironment();

    dumper->m_dumpperiod   = parser("dumpperiod").asInt();
    dumper->m_dumpdt       = parser("dumpdt").asDouble();
    dumper->m_bIO          = parser("io").asBool();
    dumper->m_bVP          = parser("vp").asBool();
    dumper->m_bHDF         = parser("hdf").asBool();
    dumper->m_bHDF_SLICE   = parser("hdf_slice").asBool();
    dumper->m_heavySkipStep= parser("heavyskipstep").asInt();
    dumper->m_channels     = parser("channels").asString();

    (*dumper)(step_id, t, NSTEPS, TEND);

    if (step_id > stepend) {
      break;
    }

    if (isroot)
      std::cout 
        << "--> t=" << t 
        << ", dt=" << dt 
        << ", step=" << step_id 
        << ", tend=" << TEND 
        << std::endl;

    dt = (*stepper)(dt, t);
    t += dt;
    ++step_id;
  }
}

#endif /* TEST_STEADYSTATEMPI_H_MT52JRKU */
