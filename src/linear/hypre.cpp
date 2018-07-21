#include "HYPRE_struct_ls.h"

#include "hypre.h"

struct Hypre::HypreData {
  HYPRE_StructGrid     grid;
  HYPRE_StructStencil  stencil;
  HYPRE_StructMatrix   a;
  HYPRE_StructVector   r;
  HYPRE_StructVector   x;
  HYPRE_StructSolver   solver;
  HYPRE_StructSolver   precond;
};

Hypre::Hypre(MPI_Comm comm, std::vector<Block> bb, MIdx gs,
             std::vector<bool> per, Scal tol, int print,
             std::string solver, int maxiter) 
  : dim(gs.size()), bb(bb), hd(new HypreData),
    solver_(solver), maxiter_(maxiter)
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
  HYPRE_StructGridCreate(comm, dim, &hd->grid);
  // Add boxes to grid
  for (auto& b : bb) {
    HYPRE_StructGridSetExtents(hd->grid, b.l.data(), b.u.data());
  }
  // Periodic
  std::vector<int> nper(per.size());
  for (size_t i = 0; i < dim; ++i) {
    nper[i] = per[i] ? gs[i] : 0;
  }
  HYPRE_StructGridSetPeriodic(hd->grid, nper.data());
  HYPRE_StructGridAssemble(hd->grid);

  // Stencil
  HYPRE_StructStencilCreate(dim, st.size(), &hd->stencil);
  for (size_t j = 0; j < st.size(); ++j) {
    HYPRE_StructStencilSetElement(hd->stencil, j, st[j].data());
  }

  // Matrix
  HYPRE_StructMatrixCreate(comm, hd->grid, hd->stencil, &hd->a);
  HYPRE_StructMatrixInitialize(hd->a);
  for (auto& b : bb) {
    std::vector<int> sti(st.size()); // stencil index (1-to-1)
    for (size_t i = 0; i < sti.size(); ++i) {
      sti[i] = (int)i;
    }

    HYPRE_StructMatrixSetBoxValues(
        hd->a, b.l.data(), b.u.data(), st.size(), sti.data(), b.a->data());
  }
  HYPRE_StructMatrixAssemble(hd->a);

  // Rhs and initial guess
  HYPRE_StructVectorCreate(comm, hd->grid, &hd->r);
  HYPRE_StructVectorCreate(comm, hd->grid, &hd->x);
  HYPRE_StructVectorInitialize(hd->r);
  HYPRE_StructVectorInitialize(hd->x);
  for (auto& b : bb) {
    HYPRE_StructVectorSetBoxValues(hd->r, b.l.data(), b.u.data(), b.r->data());
    HYPRE_StructVectorSetBoxValues(hd->x, b.l.data(), b.u.data(), b.x->data());
  }
  HYPRE_StructVectorAssemble(hd->r);
  HYPRE_StructVectorAssemble(hd->x);

  // PCG solver
  if (solver_ == "pcg") {
    HYPRE_StructPCGCreate(comm, &hd->solver);
    HYPRE_StructPCGSetMaxIter(hd->solver, maxiter_);
    HYPRE_StructPCGSetTol(hd->solver, tol); 
    //HYPRE_StructPCGSetTwoNorm(hd->solver, 1);
    //HYPRE_StructPCGSetRelChange(hd->solver, 0);
    HYPRE_StructPCGSetPrintLevel(hd->solver, print); 

    HYPRE_StructPCGSetup(hd->solver, hd->a, hd->r, hd->x);
    HYPRE_StructPCGSolve(hd->solver, hd->a, hd->r, hd->x);
  }

  // PCG solver with SMG precond
  if (solver_ == "pcg+smg") {
    // solver
    HYPRE_StructPCGCreate(comm, &hd->solver);
    HYPRE_StructPCGSetMaxIter(hd->solver, maxiter_);
    HYPRE_StructPCGSetTol(hd->solver, tol); 
    //HYPRE_StructPCGSetTwoNorm(hd->solver, 1);
    //HYPRE_StructPCGSetRelChange(hd->solver, 0);
    HYPRE_StructPCGSetPrintLevel(hd->solver, print); 

    // precond
    HYPRE_StructSMGCreate(comm, &hd->precond);
    HYPRE_StructSMGSetMemoryUse(hd->precond, 0);
    HYPRE_StructSMGSetMaxIter(hd->precond, 1);
    HYPRE_StructSMGSetTol(hd->precond, 0.0);
    HYPRE_StructSMGSetZeroGuess(hd->precond);
    HYPRE_StructSMGSetNumPreRelax(hd->precond, 1);
    HYPRE_StructSMGSetNumPostRelax(hd->precond, 1);

    // setup
    HYPRE_StructPCGSetPrecond(hd->solver, HYPRE_StructSMGSolve,
                              HYPRE_StructSMGSetup, hd->precond);
    HYPRE_StructPCGSetup(hd->solver, hd->a, hd->r, hd->x);
    HYPRE_StructPCGSolve(hd->solver, hd->a, hd->r, hd->x);
  }

  // SMG solver
  if (solver_ == "smg") {
    HYPRE_StructSMGCreate(comm, &hd->solver);
    HYPRE_StructSMGSetMemoryUse(hd->solver, 0);
    HYPRE_StructSMGSetMaxIter(hd->solver, maxiter_);
    HYPRE_StructSMGSetTol(hd->solver, tol);
    HYPRE_StructSMGSetRelChange(hd->solver, 0);
    HYPRE_StructSMGSetPrintLevel(hd->solver, print); 
    HYPRE_StructSMGSetNumPreRelax(hd->solver, 1);
    HYPRE_StructSMGSetNumPostRelax(hd->solver, 1);

    HYPRE_StructSMGSetup(hd->solver, hd->a, hd->r, hd->x);
  }

  // GMRES solver
  if (solver == "gmres") {
    HYPRE_StructGMRESCreate(comm, &hd->solver);
    HYPRE_StructGMRESSetTol(hd->solver, tol);
    HYPRE_StructGMRESSetPrintLevel(hd->solver, print);
    HYPRE_StructGMRESSetMaxIter(hd->solver, maxiter_);
    HYPRE_StructGMRESSetup(hd->solver, hd->a, hd->r, hd->x);
  }

  if (solver == "zero") {
    // nop
  }
}

void Hypre::Solve() {
  if (solver_ == "pcg+smg" || solver_ == "pcg") {
    HYPRE_StructPCGSolve(hd->solver, hd->a, hd->r, hd->x);
  }
  if (solver_ == "smg") {
    HYPRE_StructSMGSolve(hd->solver, hd->a, hd->r, hd->x);
  }
  if (solver_ == "gmres") {
    HYPRE_StructGMRESSolve(hd->solver, hd->a, hd->r, hd->x);
  }

  if (solver_ == "zero") {
    // Zero solution
    for (auto& b : bb) {
      for (auto& e : (*b.x)) {
        e = 0.;
      }
    }
  } else {
    // Copy solution
    for (auto& b : bb) {
      HYPRE_StructVectorGetBoxValues(hd->x, b.l.data(), b.u.data(), b.x->data());
    }
  }

}
Hypre::~Hypre() {
  // Destroy
  HYPRE_StructGridDestroy(hd->grid);
  HYPRE_StructStencilDestroy(hd->stencil);
  HYPRE_StructMatrixDestroy(hd->a);
  HYPRE_StructVectorDestroy(hd->r);
  HYPRE_StructVectorDestroy(hd->x);
  if (solver_ == "pcg+smg" || solver_ == "pcg") {
    HYPRE_StructPCGDestroy(hd->solver);
  }
  if (solver_ == "smg") {
    HYPRE_StructSMGDestroy(hd->solver);
  }
  if (solver_ == "gmres") {
    HYPRE_StructGMRESDestroy(hd->solver);
  }
  if (solver_ == "zero") {
    // nop
  }
}
