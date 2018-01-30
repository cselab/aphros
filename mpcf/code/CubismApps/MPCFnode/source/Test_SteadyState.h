/*
 *  Test_SteadyState.h
 *  MPCFnode
 *
 *  Created by Diego Rossinelli on 6/14/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#ifndef TEST_STEADYSTATE_H_J5GZ8VHD
#define TEST_STEADYSTATE_H_J5GZ8VHD

#include <cassert>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <ctime>

#include "Tests.h"
#include "BlockInfo.h"
#include "Types.h"
#include "OutputProcessing.h"


template <typename TGrid, typename TStepper, template <typename> class TSlice=Slice>
class Test_SteadyState: public Simulation
{
  public:

    Test_SteadyState(ArgumentParser& P) :
      grid(NULL), stepper(NULL), dumper(NULL),
      isroot(true),
      parser(P), restart_id(0), t(0.0), step_id(0)
  { }

    virtual ~Test_SteadyState()
    {
      if (this->dumper)  delete dumper;
      if (this->stepper) delete stepper;
      if (this->grid)    delete grid;
    }

    virtual void run();
    virtual void setup();

    inline void verbosity(const int v) { VERBOSITY = v; }


  protected:
    // simulation data
    TGrid * grid;

    // time stepper
    TStepper * stepper;

    // output generation
    OutputProcessing<TGrid,TSlice> * dumper;

    // class state
    bool isroot;

    // grid parameter
    int BPDX, BPDY, BPDZ;

    // stepper parameter
    int NSTEPS, SAVEPERIOD, ANALYSISPERIOD, VERBOSITY, REPORT_FREQ, REFRESHPERIOD;
    int LASTSAVE;
    Real TEND;
    bool bEXIT;

    // simulation parameter
    int restart_id;
    bool BC_PERIODIC[3];
    Real t, dt;
    int step_id;

    // helper
    ArgumentParser& parser;
    Profiler profiler;

    virtual void _setup_ic()
    {
      this->_ic();
    }

    // user info
    void _infoBoard() const;

    inline void _timestamp()
    {
      time_t rawtime;
      std::time(&rawtime);
      struct tm* timeinfo = std::localtime(&rawtime);
      char buf[256];
      std::strftime(buf, 256, "%A, %h %d %Y, %r %Z %z", timeinfo);
      std::cout << buf << std::endl;
    }


    virtual void _setEnvironment()
    {
      //parser.read_runtime_environment();

      ANALYSISPERIOD = parser("analysisperiod").asInt();
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

    virtual void _init()
    {
      if (isroot)
      {
        // print run configuration to stdout
        _infoBoard();
      }
    }
    virtual void _ic();

    virtual void _analysis()
    {
      // if (ANALYSISPERIOD != 0 && step_id%ANALYSISPERIOD == 0 && !bRESTART)
      // {
      //     profiler.push_start("ANALYSIS");
      //     profiler.pop_stop();
      // }
    }
    virtual void _pre_step() { }
    virtual void _post_step() { }
};


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SteadyState<TGrid,TStepper,TSlice>::_infoBoard() const
{
  // print all run gonfiguration settings to stdout
  // (0) my header
  std::cout << "----------------------------------------------------------------------" << std::endl;
  std::cout << "-----------------------------  INFO BOARD  ---------------------------" << std::endl;
  std::cout << "----------------------------------------------------------------------" << std::endl;
  // (1) Runtime Arguments
  std::cout << "* COMPILETIME ARGUMENTS:" << std::endl;
  std::cout << "----------------------------------------------------------------------" << std::endl;
  {
    std::cout.width(50); std::cout.fill('.');
    std::cout << std::left << "Blocksize"  << ": " << _BLOCKSIZE_ << std::endl;
    std::cout.width(50); std::cout.fill('.');
#ifdef _FLOAT_PRECISION_
    std::cout << std::left << "Float Precision"  << ": SINGLE" << std::endl;
#else
    std::cout << std::left << "Float Precision"  << ": DOUBLE" << std::endl;
#endif /* _FLOAT_PRECISION_ */
    std::cout.width(50); std::cout.fill('.');
    std::cout << std::left << "Reciprocal Precision"  << ": " << _PREC_LEVEL_ << std::endl;
    std::cout.width(50); std::cout.fill('.');
    std::cout << std::left << "Code Fusion Level"  << ": " << _MICROFUSION_ << std::endl;
    std::cout.width(50); std::cout.fill('.');
#ifdef _BGQ_
    std::cout << std::left << "Architecture"  << ": PowerPC" << std::endl;
#else
    std::cout << std::left << "Architecture"  << ": x86" << std::endl;
#endif /* _BGQ_ */
    std::cout.width(50); std::cout.fill('.');
#if defined(_QPX_) || defined(_QPXEMU_)
    std::cout << std::left << "Vectorization Support"  << ": YES" << std::endl;
#else
    std::cout << std::left << "Vectorization Support"  << ": NO" << std::endl;
#endif /* defined(_QPX_) || defined(_QPXEMU_) */
    std::cout.width(50); std::cout.fill('.');
#ifdef _USE_HDF_
    std::cout << std::left << "HDF5 Support"  << ": YES" << std::endl;
#else
    std::cout << std::left << "HDF5 Support"  << ": NO" << std::endl;
#endif /* _USE_HDF_ */
    std::cout.width(50); std::cout.fill('.');
#ifdef _WENO3_
    std::cout << std::left << "WENO Scheme"  << ": WENO3" << std::endl;
#else
    std::cout << std::left << "WENO Scheme"  << ": WENO5" << std::endl;
#endif /* _WENO3_ */
    std::cout.width(50); std::cout.fill('.');
#ifdef _NOK_
    std::cout << std::left << "K*div(u) Enabled?"  << ": NO" << std::endl;
#else
    std::cout << std::left << "K*div(u) Enabled?"  << ": YES" << std::endl;
#endif /* _NOK_ */
  }
  // (2) The Lab
  std::cout << "----------------------------------------------------------------------" << std::endl;
  std::cout << "* THE LAB:" << std::endl;
  std::cout << "----------------------------------------------------------------------" << std::endl;
  {
    Lab dummy;
    std::cout.width(50); std::cout.fill('.');
    std::cout << std::left << "Lab Name"  << ": " << dummy.name() << std::endl;
    std::cout.width(50); std::cout.fill('.');
    std::cout << std::left << "x-Periodic"  << ": " << (dummy.is_xperiodic() ? "true" : "false" ) << std::endl;
    std::cout.width(50); std::cout.fill('.');
    std::cout << std::left << "y-Periodic"  << ": " << (dummy.is_yperiodic() ? "true" : "false" ) << std::endl;
    std::cout.width(50); std::cout.fill('.');
    std::cout << std::left << "z-Periodic"  << ": " << (dummy.is_zperiodic() ? "true" : "false" ) << std::endl;
  }
  // (3) Runtime Arguments
  std::cout << "----------------------------------------------------------------------" << std::endl;
  std::cout << "* RUNTIME ARGUMENTS:" << std::endl;
  std::cout << "----------------------------------------------------------------------" << std::endl;
  parser.print_args();
  // (4) my footer
  std::cout << "----------------------------------------------------------------------" << std::endl;
}


  template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SteadyState<TGrid,TStepper,TSlice>::_setup_parameter()
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
  ANALYSISPERIOD = parser("-analysisperiod").asInt(0);
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
}


  template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SteadyState<TGrid,TStepper,TSlice>::_ic()
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
void Test_SteadyState<TGrid,TStepper,TSlice>::run()
{
  _init();

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

    (*dumper)(step_id, t, NSTEPS, TEND, profiler);

    _analysis();

    if (step_id > stepend) {
      bEXIT = true;
    }

    if (bEXIT)
    {
      break;
    }

    if (isroot && (REPORT_FREQ != 0 && step_id%REPORT_FREQ == 0))
    {
      _timestamp();
      profiler.printSummary();
    }

    if (isroot)
      std::cout 
        << "--> t=" << t 
        << ", dt=" << dt 
        << ", step=" << step_id 
        << ", tend=" << TEND 
        << std::endl;

    // (7)
    _pre_step();

    profiler.push_start("STEP");
    dt = (*stepper)(dt, t);
    profiler.pop_stop();

    t += dt;
    ++step_id;

    _post_step();

  }
}


template <typename TGrid, typename TStepper, template <typename> class TSlice>
void Test_SteadyState<TGrid,TStepper,TSlice>::setup()
{
  _setup_parameter();

  stepper = new TStepper(*grid, parser, VERBOSITY);
  assert(stepper != NULL);

  dumper = new OutputProcessing<TGrid,TSlice>(parser, *grid, isroot);
  assert(dumper != NULL);
  dumper->register_all(*grid);

  _setup_ic();
}
#endif /* TEST_STEADYSTATE_H_J5GZ8VHD */
