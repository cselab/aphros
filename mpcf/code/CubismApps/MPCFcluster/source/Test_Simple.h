#pragma once

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


template <typename TGrid, typename TStepper>
class Test_Simple : public Simulation
{
  public:
    Test_Simple(const MPI_Comm comm) :
      grid(NULL), stepper(NULL), t(0.), dt(1.),
      m_comm_world(comm)
  {
    int rank;
    MPI_Comm_rank(m_comm_world, &rank);
    isroot = (0 == rank);
  }
    virtual ~Test_Simple() {
      if (stepper) delete stepper;
      if (grid)    delete grid;
    }

    virtual void setup();


  private:
    TGrid * grid;
    TStepper * stepper;
    double dt, t;
    bool isroot;

    int BPDX, BPDY, BPDZ;

    int restart_id;
    bool BC_PERIODIC[3];

    MPI_Comm m_comm_world;

    int XPESIZE, YPESIZE, ZPESIZE;

    virtual void _setup_parameter();

    virtual void _ic();

    virtual void run();

    virtual void _init()
    {
      MPI_Barrier(m_comm_world);
    }
};


// class implementation
template <typename TGrid, typename TStepper>
void Test_Simple<TGrid,TStepper>::_setup_parameter()
{
  BPDX       = 10;
  BPDY       = 10;
  BPDZ       = 1; 
  Simulation_Environment::extent = 1.;

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

  // some checks
  assert(BPDX >= 1);
  assert(BPDY >= 1);
  assert(BPDZ >= 1);

  XPESIZE = 2;
  YPESIZE = 2;
  ZPESIZE = 1;

  const int bpdx = BPDX;
  const int bpdy = BPDY;
  const int bpdz = BPDZ;
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


template <typename TGrid, typename TStepper>
void Test_Simple<TGrid,TStepper>::setup()
{
  _setup_parameter();

  grid = new TGrid(
      XPESIZE, YPESIZE, ZPESIZE, BPDX, BPDY, BPDZ, 
      Simulation_Environment::extent, m_comm_world);

  stepper = new TStepper(*(grid));

  _ic();
}

template <typename TGrid, typename TStepper>
void Test_Simple<TGrid,TStepper>::_ic()
{
  std::vector<BlockInfo> vInfo = grid->getBlocksInfo();

  typedef typename TGrid::BlockType B;

#pragma omp parallel for
  for(int i=0; i<(int)vInfo.size(); i++)
  {
    BlockInfo info = vInfo[i];
    B& b = *(B*)info.ptrBlock;
    double* o = info.origin;

    for(int iz=0; iz<B::sizeZ; iz++)
      for(int iy=0; iy<B::sizeY; iy++)
        for(int ix=0; ix<B::sizeX; ix++)
        {
          typedef const Real CReal;
          CReal e = Simulation_Environment::extent;
          Real x = o[0] + info.h * ix / B::sizeX;
          Real y = o[1] + info.h * iy / B::sizeY;
          x /= e;
          y /= e;
          b(ix, iy, iz).alpha2 = 0.5 * (1. + std::sin(x * 10 + y * y * 5));
        }
  }
}


template <typename TGrid, typename TStepper>
void Test_Simple<TGrid,TStepper>::run()
{
  MPI_Barrier(m_comm_world);

  dt = 1.;

  for (size_t i = 0; i < 10; ++i) {
    if (isroot)
      std::cout 
        << "--> t=" << t 
        << ", dt=" << dt 
        << std::endl;

    dt = (*stepper)(dt, t);
    t += dt;
  }
}
