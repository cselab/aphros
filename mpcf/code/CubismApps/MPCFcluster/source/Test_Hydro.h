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


struct Diffusion
{
  using Lab=LabMPI;
  StencilInfo stencil;
  Real dtinvh;
  int stencil_start[3];
  int stencil_end[3];

  StencilInfo getStencil() {
    return StencilInfo(-1,-1,-1,2,2,2, true, 6, 0,1,2,3,4,6);
  }

  Diffusion(const Real _dtinvh)
    : dtinvh(_dtinvh), stencil(getStencil())
  {
    std::cerr << "Diffusion::constructor(dtinvh)" << std::endl;
    stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
    stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
  }

  Diffusion(const Diffusion& c)
    : dtinvh(c.dtinvh), stencil(c.stencil)
  {
    std::cerr << "Diffusion::copy constructor" << std::endl;
    stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
    stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
  }

  inline void operator()(Lab& lab, const BlockInfo& info, Block_t& o) const
  {
    std::cerr << "Diffusion::operator()" << std::endl;

    // At this point the only information available is:
    // - Lab, access to FluidElements via operator()(x,y,z)
    // - BlockInfo (defined in BlockInfo.h) containing index, origin, spacing
    // - Block_t=FluidBlock, 3D array with fields data and tmp
    // - dtinvh passed at construction
    //
    // Further parameters (e.g. viscosity) would be passed at construction.
    //
    // Affects: implementation of process<>()
    // Depends: interface of Kernel constructor and operator()

    // Create new instance of Kernel,
    // one instance per rank-step-block
    /*
    Kernel kernel(dtinvh);

    const Real * const srcfirst = &lab(-1,-1,-1).alpha2;
    const int labSizeRow = lab.template getActualSize<0>();
    const int labSizeSlice = labSizeRow*lab.template getActualSize<1>();
    Real * const destfirst =  &o.tmp[0][0][0][0];
    // Call kernel evaluation 
    kernel.compute(srcfirst, Block_t::gptfloats, labSizeRow, labSizeSlice,
        destfirst, Block_t::gptfloats, 
        Block_t::sizeX, Block_t::sizeX*Block_t::sizeY);
    */

  }
};


void Test_Hydro::run()
{
  MPI_Barrier(m_comm_world);

  dt = 1.;

  for (size_t i = 0; i < 10; ++i) {
    if (isroot)
      std::cerr 
        << "--> t=" << t 
        << ", dt=" << dt 
        << std::endl;

    Diffusion diffusion(dt);
    process<LabMPI>(diffusion, (GridMPI_t&)*grid, t, 0);

    //dt = (*stepper)(dt, t);
    t += dt;
  }
}
