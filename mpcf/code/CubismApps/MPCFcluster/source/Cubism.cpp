#include "Cubism.h"

#include "HYPRE_struct_ls.h"

#define NUMFR 10
#define TEND 100

using Scal = double;

Cubism::Cubism(MPI_Comm comm, KernelFactory<K>& kf, 
    int bs, Idx b, Idx p, int es, int h) 
  : bs_(bs), es_(es), h_(h), g_(p[0], p[1], p[2], b[0], b[1], b[2], 1., comm)
{
  std::vector<BlockInfo> vbi = g_.getBlocksInfo();

  #pragma omp parallel for
  for(size_t i = 0; i < vbi.size(); i++) {
    BlockInfo& bi = vbi[i];
    mk.emplace(GetIdx(bi.index), kf.Make(bi));
  }

  int r;
  MPI_Comm_rank(comm, &r);
  isroot_ = (0 == r);
}

bool Cubism::IsDone() const { 
  return step_ > TEND; 
}
void Cubism::Step() {
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

      k.ReadBuffer(l);
    }
    
    // 3. Call kernels for current stage
    for (auto& b : bb) {
      auto& k = *mk.at(GetIdx(b.index));
      k.Run();
    }

    // 4. Copy data to buffer mesh from fields collected by Comm()
    for (auto& b : bb) {
      auto& k = *mk.at(GetIdx(b.index));
      k.WriteBuffer(*(typename TGrid::BlockType*)b.ptrBlock);
    }

    // 5. Reduce
    {
      auto& f = *mk.at(GetIdx(bb[0].index)); // first kernel
      auto& vf = f.GetReduce();  // pointers to reduce

      std::vector<Scal> r(vf.size(), 0); // results

      // Check size is the same for all kernels
      for (auto& b : bb) {
        auto& k = *mk.at(GetIdx(b.index)); // kernel
        auto& v = k.GetReduce();  // pointers to reduce
        assert(v.size() == r.size());
      }
    
      // Reduce over all kernels on current rank
      for (auto& b : bb) {
        auto& k = *mk.at(GetIdx(b.index)); 
        auto& v = k.GetReduce();  
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

const std::string StreamHdf<0>::NAME = "alpha";
const std::string StreamHdf<0>::EXT = "00";

std::unique_ptr<Distr> CreateCubism(
    MPI_Comm comm, KernelFactory<Kernel>& kf, 
    int bs, Idx b, Idx p, int es, int h) {
  return unique_ptr<Distr>(new Cubism(comm, kf, bs, b, p, es, h));
}
