#pragma once

#include <cassert>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <ctime>

#include <mpi.h>

#include "BlockInfo.h"
#include "Grid.h"
#include "GridMPI.h"
#include "BlockLab.h"
#include "BlockLabMPI.h"
#include "StencilInfo.h"

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

using Real = double;

struct Elem {
  static const size_t s = 8;
  Real a[s];
  void init(Real val) { 
    for (size_t i = 0; i < s; ++i) {
      a[i] = val;
    }
  }
  void clear() {
    init(0);
  }

  Elem& operator=(const Elem&) = default;
};

struct Block {
  static const int bs = _BLOCKSIZE_;
  static const int sx = bs;
  static const int sy = bs;
  static const int sz = bs;

  // required by framework
  static const int sizeX = sx;
  static const int sizeY = sy;
  static const int sizeZ = sz;

  static const int n = sx * sy * sz;

  // floats per element
  static const int fe = sizeof(Elem) / sizeof(Real);

  using ElementType = Elem;
  using element_type = Elem;

  Elem __attribute__((__aligned__(_ALIGNBYTES_))) data[bs][bs][bs];

  Real __attribute__((__aligned__(_ALIGNBYTES_))) tmp[bs][bs][bs][fe];

  void clear_data() {
    Elem* e = &data[0][0][0];
    for(int i = 0; i < n; ++i) {
      e[i].clear();
    }
  }

  void clear_tmp() {
    Real* t = &tmp[0][0][0][0];
    for(int i = 0; i < n * fe; ++i) {
      t[i] = 0;
    }
  }

  void clear() {
    clear_data();
    clear_tmp();
  }

  inline Elem& operator()(int ix, int iy=0, int iz=0) {
    assert(ix>=0 && ix<sx);
    assert(iy>=0 && iy<sy);
    assert(iz>=0 && iz<sz);

    return data[iz][iy][ix];
  }

  /*
  template <typename Streamer>
  inline void Write(std::ofstream& output, Streamer streamer) const {
    for(int iz=0; iz<sz; iz++)
      for(int iy=0; iy<sy; iy++)
        for(int ix=0; ix<sx; ix++)
          streamer.operate(data[iz][iy][ix], output);
  }

  template <typename Streamer>
  inline void Read(std::ifstream& input, Streamer streamer) {
    for(int iz=0; iz<sz; iz++)
      for(int iy=0; iy<sy; iy++)
        for(int ix=0; ix<sx; ix++)
          streamer.operate(input, data[iz][iy][ix]);
  }
  */
};


typedef Block Block_t;  
typedef Grid<Block_t, std::allocator> GridBase;
typedef GridBase Grid_t;

template<typename BlockType, template<typename X> class Alloc=std::allocator>
class LabPer: public BlockLab<BlockType,Alloc>
{
    typedef typename BlockType::ElementType ElementTypeBlock;

public:

    virtual inline std::string name() const { return "name"; }
    bool is_xperiodic() {return true;}
    bool is_yperiodic() {return true;}
    bool is_zperiodic() {return true;}

    LabPer()
      : BlockLab<BlockType,Alloc>(){}
};

using Lab = LabPer<Block_t, std::allocator>;
typedef BlockLabMPI<Lab> LabMPI;
typedef GridMPI<Grid_t> GridMPI_t;

using TGrid = GridMPI_t;

class Test_Hydro 
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
    ~Test_Hydro() {
      if (grid)    delete grid;
    }

    void setup();
    void run();

  private:
    TGrid * grid;
    double dt, t;
    bool isroot;
    int BPDX, BPDY, BPDZ;
    int restart_id;
    bool BC_PERIODIC[3];
    MPI_Comm m_comm_world;
    int XPESIZE, YPESIZE, ZPESIZE;

    void _setup_parameter();
    void _ic();
    void _init() {
      MPI_Barrier(m_comm_world);
    }
};


// class implementation
void Test_Hydro::_setup_parameter()
{
  BPDX       = 2;
  BPDY       = 2;
  BPDZ       = 1; 

  // some post computations
  {
    Lab dummy;

    BC_PERIODIC[0] = dummy.is_xperiodic();
    BC_PERIODIC[1] = dummy.is_yperiodic();
    BC_PERIODIC[2] = dummy.is_zperiodic();
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
}


void Test_Hydro::setup()
{
  _setup_parameter();

  const Real extent = 1.;
  grid = new TGrid(
      XPESIZE, YPESIZE, ZPESIZE, BPDX, BPDY, BPDZ, 
      extent, m_comm_world);

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
          CReal e = 1.;
          Real x = o[0] + info.h * ix / B::sizeX;
          Real y = o[1] + info.h * iy / B::sizeY;
          x /= e;
          y /= e;
          b(ix, iy, iz).a[0] = 0.5 * (1. + std::sin(x * 10 + y * y * 5));
        }
  }
}


class Hydro {
 public:
   Hydro(const BlockInfo& bi) {
     name = 
         "(" + std::to_string(bi.index[0]) +
         "," + std::to_string(bi.index[1]) +
         "," + std::to_string(bi.index[2]) + ")";
   }
   void Run() {
     std::cerr << name << std::endl;
   }
   void ReadBuffer(LabMPI& l) {}
   void WriteBuffer(Block_t& o) {}
 private:
   std::string name;
};


struct Diffusion
{
  using Lab=LabMPI;
  StencilInfo stencil;
  Real dtinvh;
  using Idx = std::array<int, 3>;
  std::map<Idx, Hydro> mk;

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
    Hydro& k = i->second;
    k.Run();

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
    virtual std::unique_ptr<Kernel> Make(const BlockInfo&) = 0;
};

class HydroFactory : public KernelFactory<Hydro> {
 public:
   std::unique_ptr<Hydro> Make(const BlockInfo& bi) {
     return std::unique_ptr<Hydro>(new Hydro(bi));
   }
};

template <class Kernel>
class Distr {
 public:
  using Idx = std::array<int, 3>;

  Distr(MPI_Comm comm, KernelFactory<Kernel>& kf, 
      int bs, Idx b, Idx p, int es, int h) 
    : g_(p[0], p[1], p[2], b[0], b[1], b[2], 1., comm)
  {
    std::vector<BlockInfo> vbi = g_.getBlocksInfo();

    #pragma omp parallel for
    for(size_t i = 0; i < vbi.size(); i++)
    {
      BlockInfo& bi = vbi[i];
      mk.emplace(GetIdx(bi.index), kf.Make(bi));
    }
  }

  bool IsDone() const { return step_ > 10; }
  void Step() {
    double t = 0; // increasing timestamp
    do {
      // 1. Exchange halos in buffer mesh (by calling g_.sync())
      Diffusion rhs(g_);
      SynchronizerMPI& s = g_.sync(rhs); // only stencil needed

      LabMPI l;
      l.prepare(g_, s);

      MPI_Barrier(g_.getCartComm());

      std::vector<BlockInfo> bb = s.avail();
    
      // 2. Copy data from buffer halos to fields collected by Comm()
      for (auto& b : bb) {
        l.load(b, t);
        auto& k = mk.at(GetIdx(b.index));
        k->ReadBuffer(l);
      }
      
      // 3. Call kernels for current stage
      for (auto& b : bb) {
        auto& k = mk.at(GetIdx(b.index));
        k->Run();
      }

      // 4. Copy data to buffer mesh from fields collected by Comm()
      for (auto& b : bb) {
        auto& k = mk.at(GetIdx(b.index));
        k->WriteBuffer(*(typename TGrid::BlockType*)b.ptrBlock);
      }

      MPI_Barrier(g_.getCartComm());

      t += 1.;
      // 5. Repeat until no pending stages
    } while (t < 10.);

    ++step_;
  }

 private:
  std::map<Idx, std::unique_ptr<Kernel>> mk;

  static Idx GetIdx(const int* d) {
    return {d[0], d[1], d[2]};
  }

  TGrid g_;

  int step_ = 0;
};

void Main(MPI_Comm comm) {
  // read config files, parse arguments, maybe init global fields
  using K = Hydro;
  using KF = HydroFactory;
  using D = Distr<K>;
  using Idx = D::Idx;

  std::unique_ptr<KF> hf;

  // Kernels have to know about Distr if they want to create Stage objects
  // However, Stage can be independent on Distr.
  // Comm() should put the list of fields to exchange somewhere
  // so that Distr could do communication.
  // Comm() must be independent on implementation of Distr.
  
  Idx b{2, 2, 1}; // number of blocks 
  Idx p{2, 2, 1}; // number of ranks
  const int es = 8;
  const int h = 1;
  const int bs = 16;
  
  // Initialize buffer mesh and make Hydro for each block.
  D d(comm, *hf, bs, b, p, es, h);

  while (!d.IsDone()) {
    // At each step, first exchange halos,
    // then call Kernel::operator() for each block
    d.Step();
  }
}

/*
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
*/

#include "SynchronizerMPI.h"
#include <omp.h>
template<typename TLab, typename TKernel, typename TGrid>
inline void process(TKernel& rhs, TGrid& grid, const Real t=0.0, const bool record=false)
{
    // TKernel=Diffusion
    // TGrid=MPIGrid

    // Get synchronizer and perform communication
    SynchronizerMPI& Synch = grid.sync(rhs);

    const int nthreads = omp_get_max_threads();

    // One instance of TLab per thread
    std::vector<TLab> labs(nthreads);

    // Prepare labs: initialize boundaries and other private fields,
    // allocate memory for cacheBlock, set reference to grid.
    // CacheBlock keeps data including halos
    for(size_t i = 0; i < labs.size(); ++i)
        labs[i].prepare(grid, Synch);

    MPI_Barrier(grid.getCartComm());

    // Get list of available blocks. 
    // BlockInfo objects point to blocks in grid
    std::vector<BlockInfo> av = Synch.avail();
    const int s = av.size();
    BlockInfo * ar = &av.front();

#pragma omp parallel num_threads(nthreads)
    {
        int tid = omp_get_thread_num();
        TLab& l = labs[tid];

#pragma omp for schedule(dynamic,1)
        for(size_t i=0; i<s; i++)
        {
            // Load data from grid to lab cache (including halos)
            l.load(ar[i], t);

            // Evaluate kernel using 
            // lab cache block.data as src and
            // grid block.tmp as dst
            rhs(l, ar[i], *(typename TGrid::BlockType*)ar[i].ptrBlock);

            // The reason why blocks are copied to lab cache
            // is to have halos available.
            // BlockType normally has only inner cells.
            // Lab gives access to halo nodes with negative indices.
        }
    }

    MPI_Barrier(grid.getCartComm());
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

  MPI_Barrier(m_comm_world);
}
