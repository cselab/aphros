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
#include "BlockProcessor_MPI.h"
#include <StencilInfo.h>

#include <array>

typedef BlockLabMPI<Lab> LabMPI;
typedef GridMPI<Grid_t> GridMPI_t;
using TGrid = GridMPI_t;

class Test_Hydro : public Simulation
{
  public:
    Test_Hydro(const MPI_Comm comm) :
      grid(NULL), t(0.), dt(1.),
      m_comm_world(comm)
  {
    int rank;
    MPI_Comm_rank(m_comm_world, &rank);
    isroot = (0 == rank);
  }
    virtual ~Test_Hydro() {
      if (grid)    delete grid;
    }

    virtual void setup();


  private:
    TGrid * grid;
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
void Test_Hydro::_setup_parameter()
{
  BPDX       = 2;
  BPDY       = 2;
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


void Test_Hydro::setup()
{
  _setup_parameter();

  grid = new TGrid(
      XPESIZE, YPESIZE, ZPESIZE, BPDX, BPDY, BPDZ, 
      Simulation_Environment::extent, m_comm_world);

  // Create new instance of TStepper (e.g. HydroStep),
  // one instance per rank
  //stepper = new TStepper(*(grid));

  _ic();
}

void Test_Hydro::_ic()
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

class Kernel {
 public:
   Kernel(const BlockInfo& bi) {
     name = 
         "(" + std::to_string(bi.index[0]) +
         "," + std::to_string(bi.index[1]) +
         "," + std::to_string(bi.index[2]) + ")";
   }
   void operator()() {
     std::cerr << name << std::endl;
   }
 private:
   std::string name;
};

struct Diffusion
{
  using Lab=LabMPI;
  StencilInfo stencil;
  Real dtinvh;
  using Idx = std::array<int, 3>;
  std::map<Idx, Kernel> mk;

  static Idx GetIdx(const int* d) {
    return {d[0], d[1], d[2]};
  }

  int stencil_start[3];
  int stencil_end[3];

  StencilInfo getStencil() {
    return StencilInfo(-1,-1,-1,2,2,2, true, 6, 0,1,2,3,4,6);
  }

  Diffusion(TGrid& grid)
    : stencil(getStencil())
  {
    std::cerr << "Diffusion::constructor(dtinvh)" << std::endl;
    stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
    stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;

    std::vector<BlockInfo> vbi = grid.getBlocksInfo();

    typedef typename TGrid::BlockType B;

    #pragma omp parallel for
    for(int i=0; i<(int)vbi.size(); i++)
    {
      BlockInfo& bi = vbi[i];
      mk.emplace(GetIdx(bi.index), bi);
    }
  }

  Diffusion(const Diffusion& c) = delete;

  inline void operator()(Lab& lab, const BlockInfo& info, Block_t& o) 
  {
    if (0) {
      std::cerr 
        << "Diffusion::operator() block=(" 
        << info.index[0] << ", "
        << info.index[1] << ", "
        << info.index[2] << ")"
        << std::endl;
    }

    auto i = mk.find(GetIdx(info.index));
    assert(i != mk.end());
    Kernel& k = i->second;
    k();

    // # Current status:
    // Before each series of calls, halos are exchanged.
    // operator() is called for every block.
    // Data with halos are available in lab.
    // Results (updated values) are expected in block o.
    // Location and size of blocks is known from info.
    //
    // # Problem:
    // Persistent storage in kernels for each block.
    //
    // # Possible solution:
    // Keep a collection of objects indexed by block index
  }
};


void Test_Hydro::run()
{
  MPI_Barrier(m_comm_world);

  dt = 1.;

  Diffusion diffusion(*grid);

  for (size_t i = 0; i < 10; ++i) {
    if (isroot)
      std::cerr 
        << "--> t=" << t 
        << ", dt=" << dt 
        << std::endl;

    process<LabMPI>(diffusion, (GridMPI_t&)*grid, t, 0);

    t += dt;
  }
}
