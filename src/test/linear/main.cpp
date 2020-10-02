// Created by Petr Karnakov on 01.10.2020
// Copyright 2020 ETH Zurich

#include <iostream>

#include "debug/linear.h"
#include "distr/distrbasic.h"
#include "linear/linear.h"
#include "parse/argparse.h"
#include "solver/approx_eb.h"
#include "solver/embed.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;
using Expr = typename M::Expr;
using ExprFace = typename M::ExprFace;
using UEB = UEmbed<M>;

enum class Solver { hypre, zero, jacobi };

struct SolverInfo {
  Scal residual;
  int iter;
};

SolverInfo SolveJacobi(
    M& m, const FieldCell<Expr>& fc_system, FieldCell<Scal>& fc_sol,
    int maxiter, Scal tol) {
  auto sem = m.GetSem();
  struct {
    FieldCell<Scal> fcu;
    FieldCell<Scal> fcu_new;
    Scal maxdiff;
    int iter = 0;
    SolverInfo info;
  } * ctx(sem);
  auto& t = *ctx;
  if (sem("init")) {
    t.fcu = fc_sol;
    t.fcu_new.Reinit(m, 0);
  }
  sem.LoopBegin();
  if (sem("iter")) {
    t.maxdiff = 0;
    for (auto c : m.Cells()) {
      const auto& e = fc_system[c];
      Scal nondiag = e.back();
      for (auto q : m.Nci(c)) {
        nondiag += t.fcu[m.GetCell(c, q)] * e[1 + q];
      }
      t.fcu_new[c] = -nondiag / e[0];
      t.maxdiff = std::max(t.maxdiff, std::abs(t.fcu_new[c] - t.fcu[c]));
    }
    t.fcu.swap(t.fcu_new);
    m.Comm(&t.fcu);
    m.Reduce(&t.maxdiff, "max");
  }
  if (sem("check")) {
    t.info.residual = t.maxdiff;
    t.info.iter = t.iter;
    if (t.iter++ > maxiter || t.maxdiff < tol) {
      sem.LoopBreak();
    }
  }
  sem.LoopEnd();
  if (sem("result")) {
    fc_sol = t.fcu;
  }
  return t.info;
}

Solver GetSolver(std::string name) {
  if (name == "hypre") {
    return Solver::hypre;
  } else if (name == "zero") {
    return Solver::zero;
  } else if (name == "jacobi") {
    return Solver::jacobi;
  }
  fassert(false, "Unknown solver=" + name);
}

SolverInfo Solve(
    M& m, const FieldCell<Expr>& fc_system, FieldCell<Scal>& fc_sol,
    Solver solver, int maxiter, Scal tol) {
  switch (solver) {
    case Solver::hypre:
      Solve(fc_system, &fc_sol, fc_sol, M::LS::T::symm, m, "symm");
      return SolverInfo{m.GetResidual(), m.GetIter()};
    case Solver::zero:
      fc_sol.Reinit(m, 0);
      return SolverInfo{0, 0};
    case Solver::jacobi:
      return SolveJacobi(m, fc_system, fc_sol, maxiter, tol);
  }
  fassert(false);
}

void Run(M& m, Vars& var) {
  auto sem = m.GetSem();
  struct {
    FieldCell<Scal> fc_rhs;
    FieldCell<Scal> fc_sol;
    FieldCell<Scal> fc_sol_exact;
    FieldCell<Scal> fc_diff;
    FieldFace<Scal> ff_rho;
    FieldCell<Expr> fc_system;
    MapEmbed<BCond<Scal>> mebc;
    Solver solver;
    std::vector<generic::Vect<Scal, 3>> norms;
    SolverInfo info;
  } * ctx(sem);
  auto& t = *ctx;
  if (sem()) {
    // exact solution
    t.fc_sol_exact.Reinit(m);
    for (auto c : m.CellsM()) {
      auto x = c.center;
      t.fc_sol_exact[c] = std::sin(M_PI * std::pow(x[0], 1)) *
                          std::sin(M_PI * std::pow(x[1], 2)) *
                          std::sin(M_PI * std::pow(x[2], 3));
    }
    m.Comm(&t.fc_sol_exact);
  }
  if (sem()) {
    // resistivity
    t.ff_rho.Reinit(m);
    for (auto f : m.FacesM()) {
      t.ff_rho[f] = (Vect(0.5, 0.25, 0.125).dist(f.center()) < 0.2 ? 1000 : 1);
    }

    // system, only lhs
    FieldCell<Scal> fc_zero(m, 0);
    const auto ffg = UEmbed<M>::GradientImplicit(fc_zero, t.mebc, m);
    t.fc_system.Reinit(m, Expr::GetUnit(0));
    for (auto c : m.Cells()) {
      Expr sum(0);
      m.LoopNci(c, [&](auto q) {
        const auto cf = m.GetFace(c, q);
        const ExprFace flux = ffg[cf] / t.ff_rho[cf] * m.GetArea(cf);
        m.AppendExpr(sum, flux * m.GetOutwardFactor(c, q), q);
      });
      t.fc_system[c] = sum;
    }

    // rhs
    t.fc_rhs.Reinit(m, 0);
    for (auto c : m.Cells()) {
      t.fc_rhs[c] = UEB::Eval(t.fc_system[c], c, t.fc_sol_exact, m);
    }

    // system with rhs
    for (auto c : m.Cells()) {
      t.fc_system[c].back() = -t.fc_rhs[c];
    }

    // initial guess
    t.fc_sol.Reinit(m, 0);

    t.solver = GetSolver(var.String["solver"]);
    m.flags.linreport = 1;
  }
  if (sem.Nested("solve")) {
    t.info = Solve(
        m, t.fc_system, t.fc_sol, t.solver, var.Int["maxiter"],
        var.Double["tol"]);
  }
  if (sem("diff")) {
    t.fc_diff.Reinit(m);
    for (auto c : m.Cells()) {
      t.fc_diff[c] = t.fc_sol[c] - t.fc_sol_exact[c];
    }
  }
  if (sem.Nested("norms")) {
    t.norms = UDebug<M>::GetNorms({&t.fc_diff}, m);
  }
  if (sem("print")) {
    m.Dump(&t.fc_rhs, "rhs");
    m.Dump(&t.fc_sol, "sol");
    m.Dump(&t.fc_sol_exact, "exact");
    m.Dump(&t.fc_diff, "diff");
    if (m.IsRoot()) {
      std::cout << "\nmax_diff_exact=" << t.norms[0][2];
      std::cout << "\nresidual=" << t.info.residual;
      std::cout << "\niter=" << t.info.iter;
      std::cout << std::endl;
    }
  }
  if (sem()) {
  }
}

int main(int argc, const char** argv) {
  auto parser = ArgumentParser("Test for linear solvers.");
  parser.AddVariable<std::string>("--solver", "hypre")
      .Help("Linear solver to use. Options are: hypre, zero");
  parser.AddVariable<double>("--tol", 1e-3).Help("Convergence tolerance");
  parser.AddVariable<int>("--maxiter", 100).Help("Maximum iterations");
  auto args = parser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  std::string conf = R"EOF(
set int bx 2
set int by 2
set int bz 1

set int bsx 8
set int bsy 8
set int bsz 8

set int px 2
set int py 1
set int pz 1
)EOF";

  conf += "\nset string solver " + args.String["solver"];
  conf += "\nset double hypre_symm_tol " + args.Double.GetStr("tol");
  conf += "\nset int hypre_symm_maxiter " + args.Int.GetStr("maxiter");
  conf += "\nset double tol " + args.Double.GetStr("tol");
  conf += "\nset int maxiter " + args.Int.GetStr("maxiter");

  return RunMpiBasic<M>(argc, argv, Run, conf);
}
