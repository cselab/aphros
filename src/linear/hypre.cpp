#undef NDEBUG

#include <cassert>
#include "HYPRE_struct_ls.h"

#include "hypre.h"

struct Hypre::Imp {
  struct HypreData {
    MPI_Comm comm;
    HYPRE_StructGrid     grid;
    HYPRE_StructStencil  stencil;
    HYPRE_StructMatrix   a;
    HYPRE_StructVector   r;
    HYPRE_StructVector   x;
    HYPRE_StructSolver   solver;
    HYPRE_StructSolver   precond;
  };

  Imp(MPI_Comm comm, const std::vector<Block>& bb,
      MIdx gs, MIdx per);
  ~Imp();
  void SolverSetup(Scal tol, int print, int maxiter);
  void SolverDestroy();
  void Update();
  void Solve(Scal tol, int print, std::string solver, int maxiter);

  std::vector<Block> bb;
  HypreData hd;
  std::string solver_;
  Scal res_;
  int iter_;
};

Hypre::Imp::Imp(MPI_Comm comm, const std::vector<Block>& bb0,
                MIdx gs, MIdx per)
    : bb(bb0)
{
  assert(bb.size() > 0);
  assert(dim > 0);

  hd.comm = comm;

  // Check size of MIdx
  for (auto& b : bb) {
    assert(b.l.size() == dim);
    assert(b.u.size() == dim);
    for (auto& s : b.st) {
      if (s.size() != dim) {
        throw std::runtime_error("Hypre(): s.size() != dim");
      }
    }
  }

  // Check all blocks have the same stencil
  std::vector<MIdx> st = bb[0].st;
  for (auto& b : bb) {
    if (b.st != st) {
      throw std::runtime_error("Hypre(): b.st != st");
    }
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
  HYPRE_StructGridCreate(comm, dim, &hd.grid);
  // Add boxes to grid
  for (auto& b : bb) {
    HYPRE_StructGridSetExtents(hd.grid, b.l.data(), b.u.data());
  }
  // Periodic
  std::vector<int> nper(per.size());
  for (size_t i = 0; i < dim; ++i) {
    nper[i] = per[i] ? gs[i] : 0;
  }
  HYPRE_StructGridSetPeriodic(hd.grid, nper.data());
  HYPRE_StructGridAssemble(hd.grid);

  // Stencil
  HYPRE_StructStencilCreate(dim, st.size(), &hd.stencil);
  for (size_t j = 0; j < st.size(); ++j) {
    HYPRE_StructStencilSetElement(hd.stencil, j, st[j].data());
  }

  // Matrix
  HYPRE_StructMatrixCreate(comm, hd.grid, hd.stencil, &hd.a);
  HYPRE_StructMatrixInitialize(hd.a);
  for (auto& b : bb) {
    std::vector<int> sti(st.size()); // stencil index (1-to-1)
    for (size_t i = 0; i < sti.size(); ++i) {
      sti[i] = (int)i;
    }

    HYPRE_StructMatrixSetBoxValues(
        hd.a, b.l.data(), b.u.data(), st.size(), sti.data(), b.a->data());
  }
  HYPRE_StructMatrixAssemble(hd.a);

  // Rhs and initial guess
  HYPRE_StructVectorCreate(comm, hd.grid, &hd.r);
  HYPRE_StructVectorCreate(comm, hd.grid, &hd.x);
  HYPRE_StructVectorInitialize(hd.r);
  HYPRE_StructVectorInitialize(hd.x);
  for (auto& b : bb) {
    HYPRE_StructVectorSetBoxValues(hd.r, b.l.data(), b.u.data(), b.r->data());
    HYPRE_StructVectorSetBoxValues(hd.x, b.l.data(), b.u.data(), b.x->data());
  }
  HYPRE_StructVectorAssemble(hd.r);
  HYPRE_StructVectorAssemble(hd.x);

}

void Hypre::Imp::SolverSetup(Scal tol, int print, int maxiter) {
  // PCG solver
  if (solver_ == "pcg") {
    HYPRE_StructPCGCreate(hd.comm, &hd.solver);
    HYPRE_StructPCGSetMaxIter(hd.solver, maxiter);
    HYPRE_StructPCGSetTol(hd.solver, tol); 
    HYPRE_StructPCGSetPrintLevel(hd.solver, print); 
    HYPRE_StructPCGSetup(hd.solver, hd.a, hd.r, hd.x);
  }

  // PCG solver with SMG precond
  if (solver_ == "pcg+smg") {
    // solver
    HYPRE_StructPCGCreate(hd.comm, &hd.solver);
    HYPRE_StructPCGSetMaxIter(hd.solver, maxiter);
    HYPRE_StructPCGSetTol(hd.solver, tol); 
    HYPRE_StructPCGSetPrintLevel(hd.solver, print); 

    // precond
    HYPRE_StructSMGCreate(hd.comm, &hd.precond);
    HYPRE_StructSMGSetMemoryUse(hd.precond, 0);
    HYPRE_StructSMGSetMaxIter(hd.precond, 1);
    HYPRE_StructSMGSetTol(hd.precond, 0.0);
    HYPRE_StructSMGSetZeroGuess(hd.precond);
    HYPRE_StructSMGSetNumPreRelax(hd.precond, 1);
    HYPRE_StructSMGSetNumPostRelax(hd.precond, 1);

    // setup
    HYPRE_StructPCGSetPrecond(hd.solver, HYPRE_StructSMGSolve,
                              HYPRE_StructSMGSetup, hd.precond);
    HYPRE_StructPCGSetup(hd.solver, hd.a, hd.r, hd.x);
  }

  // SMG solver
  if (solver_ == "smg") {
    HYPRE_StructSMGCreate(hd.comm, &hd.solver);
    HYPRE_StructSMGSetMemoryUse(hd.solver, 0);
    HYPRE_StructSMGSetMaxIter(hd.solver, maxiter);
    HYPRE_StructSMGSetTol(hd.solver, tol);
    HYPRE_StructSMGSetRelChange(hd.solver, 0);
    HYPRE_StructSMGSetPrintLevel(hd.solver, print); 
    HYPRE_StructSMGSetNumPreRelax(hd.solver, 1);
    HYPRE_StructSMGSetNumPostRelax(hd.solver, 1);
    HYPRE_StructSMGSetup(hd.solver, hd.a, hd.r, hd.x);
  }

  // GMRES solver
  if (solver_ == "gmres") {
    HYPRE_StructGMRESCreate(hd.comm, &hd.solver);
    HYPRE_StructGMRESSetTol(hd.solver, tol);
    HYPRE_StructGMRESSetPrintLevel(hd.solver, print);
    HYPRE_StructGMRESSetMaxIter(hd.solver, maxiter);
    HYPRE_StructGMRESSetup(hd.solver, hd.a, hd.r, hd.x);
  }

  if (solver_ == "zero") {
    // nop
  }
}

void Hypre::Imp::SolverDestroy() {
  if (solver_ == "pcg+smg" || solver_ == "pcg") {
    HYPRE_StructPCGDestroy(hd.solver);
  }
  if (solver_ == "smg") {
    HYPRE_StructSMGDestroy(hd.solver);
  }
  if (solver_ == "gmres") {
    HYPRE_StructGMRESDestroy(hd.solver);
  }
  if (solver_ == "zero") {
    // nop
  }
}

void Hypre::Imp::Update() {
  std::vector<MIdx> st = bb[0].st;
  // Matrix
  for (auto& b : bb) {
    std::vector<int> sti(st.size()); // stencil index (1-to-1)
    for (size_t i = 0; i < sti.size(); ++i) {
      sti[i] = (int)i;
    }
    HYPRE_StructMatrixSetBoxValues(
        hd.a, b.l.data(), b.u.data(), st.size(), sti.data(), b.a->data());
  }
  HYPRE_StructMatrixAssemble(hd.a);

  // Rhs and initial guess
  for (auto& b : bb) {
    HYPRE_StructVectorSetBoxValues(hd.r, b.l.data(), b.u.data(), b.r->data());
    HYPRE_StructVectorSetBoxValues(hd.x, b.l.data(), b.u.data(), b.x->data());
  }
  HYPRE_StructVectorAssemble(hd.r);
  HYPRE_StructVectorAssemble(hd.x);
}

void Hypre::Imp::Solve(Scal tol, int print, std::string solver, int maxiter) {
  solver_ = solver;

  SolverSetup(tol, print, maxiter);

  if (solver_ == "pcg+smg" || solver_ == "pcg") {
    HYPRE_StructPCGSolve(hd.solver, hd.a, hd.r, hd.x);
    HYPRE_StructPCGGetFinalRelativeResidualNorm(hd.solver, &res_);
    HYPRE_StructPCGGetNumIterations(hd.solver, &iter_);
  }
  if (solver_ == "smg") {
    HYPRE_StructSMGSolve(hd.solver, hd.a, hd.r, hd.x);
    HYPRE_StructSMGGetFinalRelativeResidualNorm(hd.solver, &res_);
    HYPRE_StructSMGGetNumIterations(hd.solver, &iter_);
  }
  if (solver_ == "gmres") {
    HYPRE_StructGMRESSolve(hd.solver, hd.a, hd.r, hd.x);
    HYPRE_StructGMRESGetFinalRelativeResidualNorm(hd.solver, &res_);
    HYPRE_StructGMRESGetNumIterations(hd.solver, &iter_);
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
      HYPRE_StructVectorGetBoxValues(
          hd.x, b.l.data(), b.u.data(), b.x->data());
    }
  }

  SolverDestroy();
}

Hypre::Hypre(MPI_Comm comm, const std::vector<Block>& bb, MIdx gs, MIdx per)
    : imp(new Imp(comm, bb, gs, per))
{}

Hypre::~Hypre() {}

void Hypre::Update() {
  imp->Update();
}

void Hypre::Solve(Scal tol, int print, std::string solver, int maxiter) {
  imp->Solve(tol, print, solver, maxiter);
}

Hypre::Scal Hypre::GetResidual() {
  return imp->res_;
}

int Hypre::GetIter() {
  return imp->iter_;
}

Hypre::Imp::~Imp() {
  HYPRE_StructGridDestroy(hd.grid);
  HYPRE_StructStencilDestroy(hd.stencil);
  HYPRE_StructMatrixDestroy(hd.a);
  HYPRE_StructVectorDestroy(hd.r);
  HYPRE_StructVectorDestroy(hd.x);
}
