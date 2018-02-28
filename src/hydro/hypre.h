#pragma once

#include <vector>
#include <cassert>

#include <mpi.h>
#include "HYPRE_struct_ls.h"


// TODO: Convention of *_ for private variables ignored


class Hypre {
 public:
  using MIdx = std::vector<int>; 
  using Scal = double;
  struct Block { // linear system ax=r
    MIdx l; // lower corner
    MIdx u; // upper corner
    std::vector<MIdx> st; // stencil
    std::vector<Scal>* a; // matrix coeffs of size n * st.size()
    std::vector<Scal>* r; // rhs of size n
    std::vector<Scal>* x; // solution and initial guess of size n
  };

  Hypre(MPI_Comm comm, 
        std::vector<Block> bb, 
        MIdx gs /*global size*/,
        std::vector<bool> per /*periodic per direction*/,
        Scal tol /*tolerance*/,
        int print /*print level*/) 
    : dim(gs.size()), bb(bb)
  {
    assert(bb.size() > 0);
    assert(dim > 0);

    // Check size of MIdx
    for (auto& b : bb) {
      assert(b.l.size() == dim);
      assert(b.u.size() == dim);
      for (auto& s : b.st) {
        assert(s.size() == dim);
      }
    }

    // Check all blocks have the same stencil
    std::vector<MIdx> st = bb[0].st;
    for (auto& b : bb) {
      assert(b.st == st);
    }

    // Check size of data 
    for (auto& b : bb) {
      size_t n = 1;
      for (size_t i = 0; i < dim; ++i) {
        n *= b.u[i] - b.l[i] + 1;
      }
      assert(b.a->size() == n * st.size());
      assert(b.r->size() == n);
      assert(b.x->size() == n);
    }

    // Grid
    HYPRE_StructGridCreate(comm, dim, &grid);
    // Add boxes to grid
    for (auto& b : bb) {
      HYPRE_StructGridSetExtents(grid, b.l.data(), b.u.data());
    }
    // Periodic
    std::vector<int> nper(per.size());
    for (size_t i = 0; i < dim; ++i) {
      nper[i] = per[i] ? gs[i] : 0;
    }
    HYPRE_StructGridSetPeriodic(grid, nper.data());
    HYPRE_StructGridAssemble(grid);

    // Stencil
    HYPRE_StructStencilCreate(dim, st.size(), &stencil);
    for (size_t j = 0; j < st.size(); ++j) {
      HYPRE_StructStencilSetElement(stencil, j, st[j].data());
    }

    // Matrix
    HYPRE_StructMatrixCreate(comm, grid, stencil, &a);
    HYPRE_StructMatrixInitialize(a);
    for (auto& b : bb) {
      std::vector<int> sti(st.size()); // stencil index (1-to-1)
      for (int i = 0; i < sti.size(); ++i) {
        sti[i] = i;
      }

      HYPRE_StructMatrixSetBoxValues(
          a, b.l.data(), b.u.data(), st.size(), sti.data(), b.a->data());
    }
    HYPRE_StructMatrixAssemble(a);

    // Rhs and initial guess
    HYPRE_StructVectorCreate(comm, grid, &r);
    HYPRE_StructVectorCreate(comm, grid, &x);
    HYPRE_StructVectorInitialize(r);
    HYPRE_StructVectorInitialize(x);
    for (auto& b : bb) {
      HYPRE_StructVectorSetBoxValues(r, b.l.data(), b.u.data(), b.r->data());
      HYPRE_StructVectorSetBoxValues(x, b.l.data(), b.u.data(), b.x->data());
    }
    HYPRE_StructVectorAssemble(r);
    HYPRE_StructVectorAssemble(x);

    // PCG solver
    //HYPRE_StructPCGCreate(comm, &solver);
    //HYPRE_StructPCGSetTol(solver, tol); 
    //HYPRE_StructPCGSetPrintLevel(solver, print); 
    //HYPRE_StructPCGSetup(solver, a, r, x);

    // GMRES solver
    HYPRE_StructGMRESCreate(comm, &solver);
    HYPRE_StructGMRESSetTol(solver, tol);
    HYPRE_StructGMRESSetPrintLevel(solver, print);
    HYPRE_StructGMRESSetup(solver, a, r, x);
  }

  void Solve() {
    //HYPRE_StructPCGSolve(solver, a, r, x);
    HYPRE_StructGMRESSolve(solver, a, r, x);

    // Copy solution
    for (auto& b : bb) {
      HYPRE_StructVectorGetBoxValues(x, b.l.data(), b.u.data(), b.x->data());
    }
  }
  ~Hypre() {
    // Destroy
    HYPRE_StructGridDestroy(grid);
    HYPRE_StructStencilDestroy(stencil);
    HYPRE_StructMatrixDestroy(a);
    HYPRE_StructVectorDestroy(r);
    HYPRE_StructVectorDestroy(x);
    //HYPRE_StructPCGDestroy(solver);
    HYPRE_StructGMRESDestroy(solver);
  }
 
 private:
  size_t dim;
  std::vector<Block> bb;
  HYPRE_StructGrid     grid;
  HYPRE_StructStencil  stencil;
  HYPRE_StructMatrix   a;
  HYPRE_StructVector   r;
  HYPRE_StructVector   x;
  HYPRE_StructSolver   solver;
  HYPRE_StructSolver   precond;
};
