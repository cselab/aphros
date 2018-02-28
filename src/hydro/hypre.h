#pragma once

#include <vector>

#include <mpi.h>
#include "HYPRE_struct_ls.h"


class Hypre {
 public:
  using MIdx = std::vector<int>; 
  struct Block { // linear system ax=b
    MIdx l; // lower corner
    MIdx u; // upper corner
    std::vector<MIdx> st; // stencil
    std::vector<Scal>* a; // matrix coeffs of size n * st.size()
    std::vector<Scal>* b; // rhs of size n
    std::vector<Scal>* x; // solution of size n
  };

  Hypre(MPI_Comm comm, std::vector<Block> bb, MIdx gs /*global size*/) {
    auto& f = *mk.at(GetIdx(bb[0].index)); // first block
    auto& mf = f.GetMesh();
    auto& vf = mf.GetSolve();  // LS to solve

    // Check size is the same for all blocks
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

      // Grid
      HYPRE_StructGridCreate(comm, 3, &grid);

      // Add boxes to grid
      for (auto& b : bb) {
        auto d = b.index;
        using B = Block_t;
        int l[3] = {d[0] * B::sx, d[1] * B::sy, d[2] * B::sz};
        int u[3] = {l[0] + B::sx - 1, l[1] + B::sy - 1, l[2] + B::sz - 1};

        HYPRE_StructGridSetExtents(grid, l, u);
      }

      // Periodic
      if (par.Int["hypre_periodic"]) {
        int per[] = {gs[0], gs[1], gs[2]}; 
        HYPRE_StructGridSetPeriodic(grid, per);
      }

      HYPRE_StructGridAssemble(grid);

      // Stencil
      HYPRE_StructStencilCreate(3, st.size(), &stencil);
      for (size_t i = 0; i < st.size(); ++i) {
        MIdx e = st[i];
        int o[] = {e[0], e[1], e[2]};
        HYPRE_StructStencilSetElement(stencil, i, o);
      }

      // Matrix
      HYPRE_StructMatrixCreate(comm, grid, stencil, &a);
      HYPRE_StructMatrixInitialize(a);

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
        auto& m = k.GetMesh();
        auto& v = m.GetSolve();  
        auto& s = v[j]; // LS

        HYPRE_StructMatrixSetBoxValues(
            a, l, u, st.size(), sti.data(), s.a->data());
      }
      HYPRE_StructMatrixAssemble(a);

      // Rhs and initial guess
      HYPRE_StructVectorCreate(comm, grid, &b);
      HYPRE_StructVectorCreate(comm, grid, &x);

      HYPRE_StructVectorInitialize(b);
      HYPRE_StructVectorInitialize(x);

      for (auto& bi : bb) {
        auto d = bi.index;
        using B = Block_t;
        int l[3] = {d[0] * B::sx, d[1] * B::sy, d[2] * B::sz};
        int u[3] = {l[0] + B::sx - 1, l[1] + B::sy - 1, l[2] + B::sz - 1};

        auto& k = *mk.at(GetIdx(bi.index)); 
        auto& m = k.GetMesh();
        auto& v = m.GetSolve();  
        LS& s = v[j]; 

        HYPRE_StructVectorSetBoxValues(b, l, u, s.b->data());
        HYPRE_StructVectorSetBoxValues(x, l, u, s.x->data());
      }
      HYPRE_StructVectorAssemble(b);
      HYPRE_StructVectorAssemble(x);

      /*
      // PCG solver
      HYPRE_StructPCGCreate(comm, &solver);
      HYPRE_StructPCGSetTol(solver, par.Double["hypre_tol"]); 
      HYPRE_StructPCGSetPrintLevel(solver, par.Int["hypre_print"]); 
      HYPRE_StructPCGSetup(solver, a, b, x);
      HYPRE_StructPCGSolve(solver, a, b, x);
      */

      // GMRES solver
      HYPRE_StructGMRESCreate(comm, &solver);
      HYPRE_StructGMRESSetTol(solver, par.Double["hypre_tol"]);
      HYPRE_StructGMRESSetPrintLevel(solver, par.Int["hypre_print"]);
      HYPRE_StructGMRESSetup(solver, a, b, x);
      HYPRE_StructGMRESSolve(solver, a, b, x);

      // Copy solution
      for (auto& bi : bb) {
        auto d = bi.index;
        using B = Block_t;
        int l[3] = {d[0] * B::sx, d[1] * B::sy, d[2] * B::sz};
        int u[3] = {l[0] + B::sx - 1, l[1] + B::sy - 1, l[2] + B::sz - 1};

        auto& k = *mk.at(GetIdx(bi.index)); 
        auto& m = k.GetMesh();
        auto& v = m.GetSolve();  
        auto& s = v[j]; // LS

        HYPRE_StructVectorGetBoxValues(x, l, u, s.x->data());
      }

      // Destroy
      HYPRE_StructGridDestroy(grid);
      HYPRE_StructStencilDestroy(stencil);
      HYPRE_StructMatrixDestroy(a);
      HYPRE_StructVectorDestroy(b);
      HYPRE_StructVectorDestroy(x);
      //HYPRE_StructPCGDestroy(solver);
      HYPRE_StructGMRESDestroy(solver);
    }

    for (auto& b : bb) {
      auto& k = *mk.at(GetIdx(b.index)); 
      auto& m = k.GetMesh();
      m.ClearSolve();
    }
  }
  void Solve() {
  }
  ~Hypre() {}
};
