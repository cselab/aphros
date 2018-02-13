#pragma once

#include <vector>
#include <map>
#include "HYPRE_struct_ls.h"

#include "ILocal.h"
#include "../../hydro/vect.hpp"
#include "../../hydro/mesh3d.hpp"

#define NUMFR 10
#define TEND 100

using Scal = double;

template <class KF>
class Local : public Distr {
 public:
  Local(KF& kf, int bs, Idx b, Idx p, int es, int h);
  using K = typename KF::K;
  using M = typename KF::M;
  using MIdx = typename  M::MIdx;
  using Vect = typename M::Vect;
  using Rect = geom::Rect<Vect>;

  bool IsDone() const;
  void Step();

 private:
  M CreateMesh(int bs, Idx b, Idx p, int es, int h);
  M gm; // global mesh
  std::vector<MyBlockInfo> bb_;
  std::vector<geom::FieldCell<Scal>> buf_; // buffer on mesh
  std::map<Idx, std::unique_ptr<K>> mk;

  int bs_; // block size
  int es_; // element size in Scal
  int hl_; // number of halo cells (same in all directions)

  int step_ = 0;
  int stage_ = 0;
  int frame_ = 0;

  bool isroot_;
};

template <class KF>
auto Local<KF>::CreateMesh(int bs, Idx b, Idx p, int es, int hl) -> M {
  // Init global mesh
  // TODO: hl from Distr
  MIdx ms(bs_, bs_, bs_); // block size 
  MIdx mb(b[0], b[1], b[2]); // number of blocks
  MIdx mp(p[0], p[1], p[2]); // number of PEs
  MIdx mm = mp * mb * ms; // total size in cells (without halos)

  Scal ext = 1.; // TODO: extent from par
  Scal h = ext / std::max(std::max(mm[0], mm[1]), mm[2]);
  Vect d0(0); // origin coord
  Vect d1 = d0 + Vect(mm) * h;      // end coord
  Rect d(d0, d1);

  MIdx o(0); // origin index
  std::cout 
      << "o=" << o 
      << " dom=" << d0 << "," << d1 
      << " h=" << h
      << std::endl;
  
  return geom::InitUniformMesh<M>(d, o, mm);
}

template <class KF>
Local<KF>::Local(KF& kf, int bs, Idx b, Idx p, int es, int hl) 
  : bs_(bs), es_(es), hl_(hl), buf_(es_)
{
  gm = CreateMesh(bs, b, p, es, hl);

  // Resize buffer for mesh
  for (auto& u : buf_) {
    u.Reinit(gm);
  }
  
  // Fill block info
  MIdx ms(bs_, bs_, bs_); // block size 
  MIdx mb(b[0], b[1], b[2]); // number of blocks
  MIdx mp(p[0], p[1], p[2]); // number of PEs
  geom::BlockCells<3> bc(mb * mp);
  using geom::IdxNode;
  Scal h = gm.GetNode(IdxNode(1)) - gm.GetNode(IdxNode(0));
  std::cerr << "h from gm = " << h << std::endl;
  for (MIdx i : bc) {
    MyBlockInfo b;
    IdxNode n = gm.GetBlockNodes().GetIdx(i * ms);
    Vect o = gm.GetCenter(n);
    std::cerr << "o=" << o << " n=" << n <<  " i=" << i << std::endl;
    for (int q = 0; q < 3; ++q) {
      b.index[q] = i[q];
      b.origin[q] = i[q];
    }
    b.h_gridpoint = h;
    b.ptrBlock = nullptr;
    bb_.push_back(b);
  }

  for(size_t i = 0; i < bb_.size(); i++) {
    MyBlockInfo& bi = bb_[i];
    auto up = kf.Make(bi);
    mk.emplace(GetIdx(bi.index), std::unique_ptr<K>(dynamic_cast<K*>(up.release())));
  }
}

template <class KF>
bool Local<KF>::IsDone() const { 
  return step_ > TEND; 
}

template <class KF>
void Local<KF>::Step() {
  if (isroot_) {
    std::cerr << "***** STEP " << step_ << " ******" << std::endl;
  }
  do {
    if (isroot_) {
      std::cerr << "*** STAGE abs=" << stage_ << " ***" << std::endl;
    }
    // 1. Exchange halos in buffer mesh
    // Do communication and get all blocks
    std::vector<MyBlockInfo> bb; // TODO: fill
    //assert(!bb.empty());
  
    // 2. Copy data from buffer halos to fields collected by Comm()
    for (auto& b : bb) {
      auto& k = *mk.at(GetIdx(b.index)); // kernel
      auto& m = k.GetMesh();
      //ReadBuffer(m);
    }
    
    // 3. Call kernels for current stage
    for (auto& b : bb) {
      auto& k = *mk.at(GetIdx(b.index));
      k.Run();
    }

    // 4. Copy data to buffer mesh from fields collected by Comm()
    for (auto& b : bb) {
      auto& k = *mk.at(GetIdx(b.index)); // kernel
      auto& m = k.GetMesh();
      //WriteBuffer(m);
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

      // Write results to all kernels on current rank
      for (auto& b : bb) {
        auto& k = *mk.at(GetIdx(b.index)); 
        auto& m = k.GetMesh();
        auto& v = m.GetReduce();  
        for (size_t i = 0; i < r.size(); ++i) {
          *v[i] = r[i];
        }
      }

      // Clear reduce requests
      for (auto& b : bb) {
        auto& k = *mk.at(GetIdx(b.index)); 
        auto& m = k.GetMesh();
        m.ClearReduce();
      }
    }

    // 6. Solve 
    {
      MPI_Comm comm = MPI_COMM_WORLD;
      auto& f = *mk.at(GetIdx(bb[0].index)); // first kernel
      auto& mf = f.GetMesh();
      auto& vf = mf.GetSolve();  // LS to solve

      // Check size is the same for all kernels
      for (auto& b : bb) {
        auto& k = *mk.at(GetIdx(b.index)); // kernel
        auto& m = k.GetMesh();
        auto& v = m.GetSolve();  // pointers to reduce
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
          using B = MyBlock;
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
          using B = MyBlock;
          int l[3] = {d[0] * B::sx, d[1] * B::sy, d[2] * B::sz};
          int u[3] = {l[0] + B::sx - 1, l[1] + B::sy - 1, l[2] + B::sz - 1};

          std::vector<int> sti(st.size()); // stencil index (1-1)
          for (int i = 0; i < sti.size(); ++i) {
            sti[i] = i;
          }

          auto& k = *mk.at(GetIdx(b.index)); 
          auto& m = k.GetMesh();
          auto& v = m.GetSolve();  
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
          using B = MyBlock;
          int l[3] = {d[0] * B::sx, d[1] * B::sy, d[2] * B::sz};
          int u[3] = {l[0] + B::sx - 1, l[1] + B::sy - 1, l[2] + B::sz - 1};

          auto& k = *mk.at(GetIdx(bi.index)); 
          auto& m = k.GetMesh();
          auto& v = m.GetSolve();  
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
          using B = MyBlock;
          int l[3] = {d[0] * B::sx, d[1] * B::sy, d[2] * B::sz};
          int u[3] = {l[0] + B::sx - 1, l[1] + B::sy - 1, l[2] + B::sz - 1};

          auto& k = *mk.at(GetIdx(bi.index)); 
          auto& m = k.GetMesh();
          auto& v = m.GetSolve();  
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

      for (auto& b : bb) {
        auto& k = *mk.at(GetIdx(b.index)); 
        auto& m = k.GetMesh();
        m.ClearSolve();
      }
    }


    stage_ += 1;

    // 6. Check for pending stages
    {
      int np = 0;
      for (auto& b : bb) {
        auto& k = *mk.at(GetIdx(b.index));
        auto& m = k.GetMesh();
        if (m.Pending()) {
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
    //DumpHDF5_MPI<TGrid, StreamHdf<0>>(g_, frame_, step_*1., "p" + suff);
    ++frame_;
  }
  ++step_;
}

