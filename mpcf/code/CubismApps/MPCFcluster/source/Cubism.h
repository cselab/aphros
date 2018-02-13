#pragma once

#include "BlockInfo.h"
#include "Grid.h"
#include "GridMPI.h"
#include "BlockLab.h"
#include "BlockLabMPI.h"
#include "StencilInfo.h"
#include "HDF5Dumper_MPI.h"
#include "ICubism.h"

#include "HYPRE_struct_ls.h"

#define NUMFR 10
#define TEND 100

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
  static const int n = sx * sy * sz;

  // required by framework
  static const int sizeX = sx;
  static const int sizeY = sy;
  static const int sizeZ = sz;


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

using Scal = double;

template <class KF>
class Cubism {
 public:
  using Idx = std::array<int, 3>;

  Cubism(MPI_Comm comm, KF& kf, 
      int bs, Idx b, Idx p, int es, int h);
  using K = typename KF::K;

  bool IsDone() const;
  void Step();

 private:
  std::map<Idx, std::unique_ptr<K>> mk;

  static Idx GetIdx(const int* d) {
    return {d[0], d[1], d[2]};
  }

  int bs_; // block size
  int es_; // element size in Scal
  int h_; // number of halo cells (same in all directions)

  TGrid g_;

  int step_ = 0;
  int stage_ = 0;
  int frame_ = 0;

  static StencilInfo GetStencil(int h) {
    return StencilInfo(-h,-h,-h,h+1,h+1,h+1, true, 8, 0,1,2,3,4,5,6,7);
  }

  bool isroot_;

  //void ReadBuffer(Hydro<MeshStructured>&);
};


template <int ID>
struct StreamHdf {
  static const std::string NAME;
  static const std::string EXT;
  static const int NCHANNELS = 1;
  static const int CLASS = 0;
  struct T { Real a[8]; };

  using B = Block;
  B& b;

  StreamHdf(B& b): b(b) {}

  // write
  void operate(const int ix, const int iy, const int iz, Real out[0]) const
  {
    const T& in = *((T*)&b.data[iz][iy][ix].a[0]);
    out[0] = in.a[ID];
  }

  // read
  void operate(const Real out[0], const int ix, const int iy, const int iz) const
  {
    T& in = *((T*)&b.data[iz][iy][ix].a[0]);
    in.a[ID] = out[0];
  }

  static const char * getAttributeName() { return "Scalar"; }
};

struct FakeProc {
  StencilInfo stencil;
  explicit FakeProc(StencilInfo si) 
    : stencil(si)
  {}
};

template <class KF>
Cubism<KF>::Cubism(MPI_Comm comm, KF& kf, 
    int bs, Idx b, Idx p, int es, int h) 
  : bs_(bs), es_(es), h_(h), g_(p[0], p[1], p[2], b[0], b[1], b[2], 1., comm)
{
  std::vector<BlockInfo> vbi = g_.getBlocksInfo();

  #pragma omp parallel for
  for(size_t i = 0; i < vbi.size(); i++) {
    BlockInfo& bi = vbi[i];
    MyBlockInfo mbi;
    for (int j = 0; j < 3; ++j) {
      mbi.index[j] = bi.index[j];
      mbi.origin[j] = bi.origin[j];
    }
    mbi.h_gridpoint = bi.h_gridpoint;
    mbi.ptrBlock = bi.ptrBlock;
    mk.emplace(GetIdx(bi.index), kf.Make(mbi));
  }

  int r;
  MPI_Comm_rank(comm, &r);
  isroot_ = (0 == r);
}

template <class KF>
bool Cubism<KF>::IsDone() const { 
  return step_ > TEND; 
}

template <class KF>
void Cubism<KF>::Step() {
  MPI_Barrier(g_.getCartComm());
  if (isroot_) {
    std::cerr << "***** STEP " << step_ << " ******" << std::endl;
  }
  do {
    MPI_Barrier(g_.getCartComm());
    if (isroot_) {
      std::cerr << "*** STAGE abs=" << stage_ << " ***" << std::endl;
    }
    // 1. Exchange halos in buffer mesh
    FakeProc fp(GetStencil(h_));       // object with field 'stencil'
    SynchronizerMPI& s = g_.sync(fp); 

    LabMPI l;
    l.prepare(g_, s);   // allocate memory for lab cache

    MPI_Barrier(g_.getCartComm());

    // Do communication and get all blocks
    std::vector<BlockInfo> bb = s.avail();
    assert(!bb.empty());
  
    // 2. Copy data from buffer halos to fields collected by Comm()
    for (auto& b : bb) {
      l.load(b, stage_);
      MPI_Barrier(g_.getCartComm());
      auto& k = *mk.at(GetIdx(b.index));

      //k.ReadBuffer(l);
    }
    
    // 3. Call kernels for current stage
    for (auto& b : bb) {
      auto& k = *mk.at(GetIdx(b.index));
      k.Run();
    }

    // 4. Copy data to buffer mesh from fields collected by Comm()
    for (auto& b : bb) {
      auto& k = *mk.at(GetIdx(b.index));
      //k.WriteBuffer(*(typename TGrid::BlockType*)b.ptrBlock);
    }

    // 5. Reduce
    {
      auto& f = *mk.at(GetIdx(bb[0].index)); // first kernel
      auto& mf = f.GetMesh();
      auto& vf = mf.GetReduce();  // pointers to reduce

      std::vector<Scal> r(vf.size(), 0); // results

      // Check size is the same for all kernels
      for (auto& b : bb) {
        auto& k = *mk.at(GetIdx(b.index)); // kernel
        auto& m = k.GetMesh();
        auto& v = m.GetReduce();  // pointers to reduce
        assert(v.size() == r.size());
      }
    
      // Reduce over all kernels on current rank
      for (auto& b : bb) {
        auto& k = *mk.at(GetIdx(b.index)); 
        auto& m = k.GetMesh();
        auto& v = m.GetReduce();  
        for (size_t i = 0; i < r.size(); ++i) {
          r[i] += *v[i];
        }
      }

      // Reduce over all ranks
      MPI_Allreduce(
          MPI_IN_PLACE, r.data(), r.size(), 
          MPI_DOUBLE, MPI_SUM, g_.getCartComm()); // TODO: type from Scal

      // Write results to all kernels on current rank
      for (auto& b : bb) {
        auto& k = *mk.at(GetIdx(b.index)); 
        auto& v = k.GetReduce();  
        for (size_t i = 0; i < r.size(); ++i) {
          *v[i] = r[i];
        }
      }
    }

    // 6. Solve 
    {
      MPI_Comm comm = g_.getCartComm();
      auto& f = *mk.at(GetIdx(bb[0].index)); // first kernel
      auto& vf = f.GetSolve();  // LS to solve

      // Check size is the same for all kernels
      for (auto& b : bb) {
        auto& k = *mk.at(GetIdx(b.index)); // kernel
        auto& v = k.GetSolve();  // pointers to reduce
        assert(v.size() == vf.size());
      }

      for (size_t j = 0; j < vf.size(); ++j) {
        using MIdx = typename K::MIdx;
        std::vector<MIdx> st = vf[j].st; // stencil

        HYPRE_StructGrid     grid;
        HYPRE_StructStencil  stencil;
        HYPRE_StructMatrix   a;
        HYPRE_StructVector   b;
        HYPRE_StructVector   x;
        HYPRE_StructSolver   solver;
        HYPRE_StructSolver   precond;

        // Create empty 3D grid object
        HYPRE_StructGridCreate(comm, 3, &grid);

        // Add boxes to grid
        for (auto& b : bb) {
          auto d = b.index;
          using B = Block_t;
          int l[3] = {d[0] * B::sx, d[1] * B::sy, d[2] * B::sz};
          int u[3] = {l[0] + B::sx - 1, l[1] + B::sy - 1, l[2] + B::sz - 1};

          HYPRE_StructGridSetExtents(grid, l, u);
        }

        // Assemble grid
        HYPRE_StructGridAssemble(grid);

        // Create emtpty 3D stencil
        HYPRE_StructStencilCreate(3, st.size(), &stencil);

        // Assign stencil entries 
        for (size_t i = 0; i < st.size(); ++i) {
          MIdx e = st[i];
          int o[] = {e[0], e[1], e[2]};
          HYPRE_StructStencilSetElement(stencil, i, o);
        }

        // Create empty matrix object 
        HYPRE_StructMatrixCreate(comm, grid, stencil, &a);

        // Indicate that the matrix coefficients are ready to be set 
        HYPRE_StructMatrixInitialize(a);

        // Set matrix coefficients
        for (auto& b : bb) {
          auto d = b.index;
          using B = Block_t;
          int l[3] = {d[0] * B::sx, d[1] * B::sy, d[2] * B::sz};
          int u[3] = {l[0] + B::sx - 1, l[1] + B::sy - 1, l[2] + B::sz - 1};

          std::vector<int> sti(st.size()); // stencil index (1-1)
          for (int i = 0; i < sti.size(); ++i) {
            sti[i] = i;
          }

          auto& k = *mk.at(GetIdx(b.index)); 
          auto& v = k.GetSolve();  
          auto& s = v[j]; // LS

          HYPRE_StructMatrixSetBoxValues(
              a, l, u, st.size(), sti.data(), s.a->data());
        }

        // Assemble matrix
        HYPRE_StructMatrixAssemble(a);

        /* Create an empty vector object */
        HYPRE_StructVectorCreate(comm, grid, &b);
        HYPRE_StructVectorCreate(comm, grid, &x);

        /* Indicate that the vector coefficients are ready to be set */
        HYPRE_StructVectorInitialize(b);
        HYPRE_StructVectorInitialize(x);

        /* Set the vector coefficients over the gridpoints in my first box */
        for (auto& bi : bb) {
          auto d = bi.index;
          using B = Block_t;
          int l[3] = {d[0] * B::sx, d[1] * B::sy, d[2] * B::sz};
          int u[3] = {l[0] + B::sx - 1, l[1] + B::sy - 1, l[2] + B::sz - 1};

          auto& k = *mk.at(GetIdx(bi.index)); 
          auto& v = k.GetSolve();  
          auto& s = v[j]; // LS

          HYPRE_StructVectorSetBoxValues(b, l, u, s.b->data());
          HYPRE_StructVectorSetBoxValues(x, l, u, s.x->data());
       }

        /* This is a collective call finalizing the vector assembly.
           The vectors are now ``ready to be used'' */
        HYPRE_StructVectorAssemble(b);
        HYPRE_StructVectorAssemble(x);


        /*
     // 6.5. Set up and use a solver (See the Reference Manual for descriptions
        //of all of the options.)
        // Create an empty PCG Struct solver 
        HYPRE_StructPCGCreate(comm, &solver);

        // Set PCG parameters
        HYPRE_StructPCGSetTol(solver, 1.0e-06);
        HYPRE_StructPCGSetPrintLevel(solver, 2);
        HYPRE_StructPCGSetMaxIter(solver, 50);

        // Use symmetric SMG as preconditioner 
        HYPRE_StructSMGCreate(comm, &precond);
        HYPRE_StructSMGSetMaxIter(precond, 1);
        HYPRE_StructSMGSetTol(precond, 0.0);
        HYPRE_StructSMGSetZeroGuess(precond);
        HYPRE_StructSMGSetNumPreRelax(precond, 1);
        HYPRE_StructSMGSetNumPostRelax(precond, 1);

        // Set preconditioner and solve
        HYPRE_StructPCGSetPrecond(solver, HYPRE_StructSMGSolve,
                                  HYPRE_StructSMGSetup, precond);
        HYPRE_StructPCGSetup(solver, a, b, x);
        HYPRE_StructPCGSolve(solver, a, b, x);

        */

    /* Create an empty PCG Struct solver */
    HYPRE_StructPCGCreate(comm, &solver);

    /* Set some parameters */
    HYPRE_StructPCGSetTol(solver, 1.0e-06); /* convergence tolerance */
    HYPRE_StructPCGSetPrintLevel(solver, 0); /* amount of info. printed */

    /* Setup and solve */
    HYPRE_StructPCGSetup(solver, a, b, x);
    HYPRE_StructPCGSolve(solver, a, b, x);

        for (auto& bi : bb) {
          auto d = bi.index;
          using B = Block_t;
          int l[3] = {d[0] * B::sx, d[1] * B::sy, d[2] * B::sz};
          int u[3] = {l[0] + B::sx - 1, l[1] + B::sy - 1, l[2] + B::sz - 1};

          auto& k = *mk.at(GetIdx(bi.index)); 
          auto& v = k.GetSolve();  
          auto& s = v[j]; // LS

          HYPRE_StructVectorGetBoxValues(x, l, u, s.x->data());
       }


     HYPRE_StructGridDestroy(grid);
     HYPRE_StructStencilDestroy(stencil);
     HYPRE_StructMatrixDestroy(a);
     HYPRE_StructVectorDestroy(b);
     HYPRE_StructVectorDestroy(x);
     HYPRE_StructPCGDestroy(solver);


     }

    }


    MPI_Barrier(g_.getCartComm());

    stage_ += 1;

    // 6. Check for pending stages
    {
      int np = 0;
      for (auto& b : bb) {
        auto& k = *mk.at(GetIdx(b.index));
        if (k.Pending()) {
          ++np;
        }
      }
      // Check either all done or all pending
      assert(np == 0 || np == bb.size());

      // Break if no pending stages
      if (!np) {
        break;
      }
    }
  } while (true);

  if (step_ % (TEND / NUMFR)  == 0) {
    auto suff = "_" + std::to_string(frame_);
    DumpHDF5_MPI<TGrid, StreamHdf<0>>(g_, frame_, step_*1., "p" + suff);
    ++frame_;
  }
  ++step_;
}

template <int i>
const std::string StreamHdf<i>::NAME = "alpha";
template <int i>
const std::string StreamHdf<i>::EXT = "00";
