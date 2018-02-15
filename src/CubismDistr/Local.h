#pragma once

#include <vector>
#include <map>
#include <mpi.h>
#include "HYPRE_struct_ls.h"

#include "ILocal.h"
#include "hydro/vect.hpp"
#include "hydro/mesh3d.hpp"
#include "hydro/output.hpp"
#include "hydro/output_paraview.hpp"

#define NUMFR 10
#define TEND 100

using Scal = double;

template <class KF>
class Local : public Distr {
 public:
  Local(MPI_Comm comm, KF& kf, int bs, Idx b, Idx p, int es, int h);
  using K = typename KF::K;
  using M = typename KF::M;
  using MIdx = typename  M::MIdx;
  using Vect = typename M::Vect;
  using Rect = geom::Rect<Vect>;
  using IdxCell = geom::IdxCell;

  bool IsDone() const;
  void Step();

 private:
  MPI_Comm comm_;
  int bs_; // block size
  int es_; // element size in Scal
  int hl_; // number of halo cells (same in all directions)
  std::vector<geom::FieldCell<Scal>> buf_; // buffer on mesh

  M gm; // global mesh
  std::unique_ptr<output::Session> session_;
  std::vector<MyBlockInfo> bb_;
  std::map<Idx, std::unique_ptr<K>> mk;

  void ReadBuffer(M& m);
  void WriteBuffer(M& m);
  M CreateMesh(int bs, Idx b, Idx p, int es, int h);

  int step_ = 0;
  int stage_ = 0;
  int frame_ = 0;

  bool isroot_ = true;
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
  
  return geom::InitUniformMesh<M>(d, o, mm, 0);
}

template <class KF>
Local<KF>::Local(MPI_Comm comm, KF& kf, int bs, Idx b, Idx p, int es, int hl) 
  : comm_(comm), bs_(bs), es_(es), hl_(hl), buf_(es_)
{
  gm = CreateMesh(bs, b, p, es, hl);

  output::Content content = {
      std::make_shared<output::EntryFunction<Scal, IdxCell, M>>(
          "p", gm, [this](IdxCell i) { return buf_[0][i]; })
      };

  session_.reset(new output::SessionParaviewStructured<M>(
          content, "title", "p", gm));

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
  Scal h = (gm.GetNode(IdxNode(1)) - gm.GetNode(IdxNode(0)))[0];
  assert(h > 0);
  std::cerr << "h from gm = " << h << std::endl;
  for (MIdx i : bc) {
    MyBlockInfo b;
    IdxNode n = gm.GetBlockNodes().GetIdx(i * ms);
    Vect o = gm.GetNode(n);
    std::cerr << "o=" << o << " n=" << n.GetRaw() <<  " i=" << i << std::endl;
    for (int q = 0; q < 3; ++q) {
      b.index[q] = i[q];
      b.origin[q] = o[q];
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
    auto& bb = bb_; // TODO: fill
    assert(!bb.empty());
  
    // 2. Copy data from buffer halos to fields collected by Comm()
    for (auto& b : bb) {
      auto& k = *mk.at(GetIdx(b.index)); // kernel
      auto& m = k.GetMesh();
      ReadBuffer(m);
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
      WriteBuffer(m);
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
      MPI_Comm comm = comm_;
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
    //auto suff = "_" + std::to_string(frame_);
    std::cerr << "Output" << std::endl;
    session_->Write(step_*1., "title:0");
    ++frame_;
  }
  ++step_;
}

template <class KF>
void Local<KF>::ReadBuffer(M& m) {
  int e = 0; // buffer field idx

  for (auto u : m.GetComm()) {
    for (auto i : m.AllCells()) {
      auto& bc = m.GetBlockCells();
      auto& gbc = gm.GetBlockCells();
      MIdx gs = gbc.GetDimensions();
      auto d = bc.GetMIdx(i);
      // periodic
      for (int j = 0; j < 3; ++j) {
        d[j] = (d[j] + gs[j]) % gs[j];
      }
      auto gi = gbc.GetIdx(d);
      (*u)[i] = buf_[e][gi];
    }
    ++e;
  }

  m.ClearComm();
}

template <class KF>
void Local<KF>::WriteBuffer(M& m) {
  using MIdx = typename M::MIdx;

  // Check buffer has enough space for all fields
  assert(m.GetComm().size() <= buf_.size() && "Too many fields for Comm()");

  int e = 0; // buffer field idx

  for (auto u : m.GetComm()) {
    for (auto i : m.Cells()) {
      auto& bc = m.GetBlockCells();
      auto& gbc = gm.GetBlockCells();
      auto d = bc.GetMIdx(i); 
      auto gi = gbc.GetIdx(d); 
      buf_[e][gi] = (*u)[i];
    }
    ++e;
  }
}
