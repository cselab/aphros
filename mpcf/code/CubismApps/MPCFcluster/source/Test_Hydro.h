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

/*
 
Basic Rules:

1. Program at interface.
First describe the interface in a separate class, then implement.

2. Use templates when:
- 

3. Use inheritance when:
- 
 
*/

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
    //
    // # Outcome:
    // Diffusion can instantiate Kernel passing BlockInfo to it.
    // 
    // Workflow:
    // Pass mesh size, block size, stencil to Cubism (here Diffusion).
    // Cubism initializes the mesh with generic elements (buffers)
    // and instantiates Kernel for each block.
    // At each simulation step (series of calls) 
    // it exhanges halos and calls Kernel::operator() of every block.
    // TODO: Kernel::operator() needs arguments for data (Lab and Block_t)
    //
    // How to pass other arguments to Kernel (e.g. P_double)?
  }
};

// Plan
// 1. Implement Comm() under assumption of collective requests
//

// Should I use a factory to instantiate Kernels?
//   Workflow:
// - Constructor of Kernel is called once per rank
// - That constructor initializes the environment 
//   reading files and parsing arguments if needed
// - Distr never calls Kernel constructor directly.
//   Instead, it assumes that Kernel has method Make(BlockInfo)
//   which creates a clone of that kernel but for a different block
// - Distr is created for that factory.
//   That's the point. Distr is a template with argument Kernel.
//   But it never instantiates Kernel, instead,
//   it calls the factory.
//
// Comm() puts fields in a list separately for each block.
// Kernel has methods WriteBuffer(), ReadBuffer()
// which write to buffer mesh after kernel call
// and read from buffer before kernel call.

template <class Kernel>
class KernelFactory {
  public:
    virtual std::unique_ptr<Kernel> Make() = 0;
};

template <class Kernel>
class Distr {
 public:
  using MIdx = std::array<size_t, 3>;
  using Grid = GridMPI<BaseGrid>;
  using Idx = std::array<int, 3>;

  Distr(const KernelFactory<Kernel>& kf, 
      size_t bs, MIdx b, MIdx p, size_t es, size_t h) 
    : g_(bs, b, p, es, h)
  {
    std::vector<BlockInfo> vbi = g_.getBlocksInfo();

    #pragma omp parallel for
    for(size_t i = 0; i < vbi.size(); i++)
    {
      BlockInfo& bi = vbi[i];
      mk.emplace(GetIdx(bi.index), kf.Make(bi));
    }
  }

  bool IsDone() const { return false; }
  void Step() {
    do {
      // 1. Exchange halos in buffer mesh (by calling g_.sync())
      bb = g_.sync().avail();
    
      // 2. Copy data from buffer halos to fields collected previously by Comm()
      for (auto& b : bb) {
        b.ReadBuffer(lab);
      }
      
      // 3. Call kernels for current stage
      for (auto& b : bb) {
        b.Run(lab);
      }

      // 4. Copy data to buffer mesh for fields collected by Comm()
      for (auto& b : bb) {
        b.WriteBuffer(b.ptr);
      }

      // 5. Repeat until no pending stages
    } while (true);
  }

 private:
  std::map<Idx, Kernel> mk;

  static Idx GetIdx(const int* d) {
    return {d[0], d[1], d[2]};
  }

  Grid g_;
};

template <class T>
void Main(int argc, char** argv) {
  using HydroFactory = T;
  using Hydro = T;
  // read config files, parse arguments, maybe init global fields
  HydroFactory hf(argc, argv);

  // This would create a kernel for a single block
  // i.e. create mesh, initialize local fields.
  // But normally it is called from a Distr.
  //hf.Make(BlockInfo(...));

  // Kernels have to know about Distr if they want to create Stage objects
  // However, Stage can be independent on Distr.
  // Comm() should put the list of fields to exchange somewhere
  // so that Distr could do communication.
  // Comm() must be independent on implementation of Distr.
  
  MIdx b(10, 10, 10); // number of blocks 
  
  // Initialize buffer mesh and make Hydro for each block.
  Distr<Hydro> d(hf, b);

  while (!d.IsDone()) {
    // At each step, first exchange halos,
    // then call Kernel::operator() for each block
    d.Step();
  }
}

template <class T>
void Example() {
  using Cubism = typename T::Cubism;
  // Initialize mesh:
  // bs : block size
  // bx, by, bz : number of blocks per PE
  // px, py, pz : number of PEs
  // es : element size
  // hl : number of halo cells (same in all dimensions)
  // Instantiate Kernel for each block.
  Cubism<Kernel> c(bs, bx, by, bz, px, py, pz, es, hl);
  while (!c.IsDone()) {
    // At each step, first exchange halos,
    // then call Kernel::operator() for each block
    c.Step();
  }
}

template <class T>
class Kern() {
  Distr d;
  std::vector<double> a;
  void Step() {
    st = d.GetStage();
    if (st()) {
      st.Comm(a);
    }
  }
}


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
