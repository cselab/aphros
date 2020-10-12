// Created by Petr Karnakov on 01.10.2020
// Copyright 2020 ETH Zurich

#include <iostream>

#include "debug/linear.h"
#include "distr/distrbasic.h"
#include "linear/hypre.h"
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

enum class Solver { hypre, zero, jacobi, conjugate };

struct SolverInfo {
  Scal residual;
  int iter;
};

template <class M>
SolverInfo SolveHypre(
    const FieldCell<typename M::Expr>& fc_system,
    const FieldCell<typename M::Scal>* fc_init,
    FieldCell<typename M::Scal>& fc_sol, M& m, int maxiter, Scal tol) {
  auto sem = m.GetSem(__func__);
  struct {
    std::unique_ptr<linear::Solver<M>> solver;
  } * ctx(sem);
  auto& t = *ctx;
  if (sem("init")) {
    typename linear::Solver<M>::Conf conf;
    typename linear::SolverHypre<M>::Extra extra;

    conf.tol = tol;
    conf.maxiter = maxiter;
    extra.solver = "pcg";

    t.solver = std::make_unique<linear::SolverHypre<M>>(conf, extra);
  }
  if (sem.Nested("solve")) {
    auto info = t.solver->Solve(fc_system, fc_init, fc_sol, m);
    return {info.residual, info.iter};
  }
  return {};
}

SolverInfo SolveConjugate(
    M& m, const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
    FieldCell<Scal>& fc_sol, int maxiter, Scal tol) {
  auto sem = m.GetSem(__func__);
  struct {
    std::unique_ptr<linear::Solver<M>> solver;
  } * ctx(sem);
  auto& t = *ctx;
  if (sem("init")) {
    typename linear::Solver<M>::Conf conf;
    typename linear::SolverConjugate<M>::Extra extra;

    conf.tol = tol;
    conf.maxiter = maxiter;

    t.solver = std::make_unique<linear::SolverConjugate<M>>(conf, extra);
  }
  if (sem.Nested("solve")) {
    auto info = t.solver->Solve(fc_system, fc_init, fc_sol, m);
    return {info.residual, info.iter};
  }
  return {};
}

SolverInfo SolveJacobi(
    M& m, const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
    FieldCell<Scal>& fc_sol, int maxiter, Scal tol) {
  auto sem = m.GetSem(__func__);
  struct {
    std::unique_ptr<linear::Solver<M>> solver;
  } * ctx(sem);
  auto& t = *ctx;
  if (sem("init")) {
    typename linear::Solver<M>::Conf conf;
    typename linear::SolverJacobi<M>::Extra extra;

    conf.tol = tol;
    conf.maxiter = maxiter;

    t.solver = std::make_unique<linear::SolverJacobi<M>>(conf, extra);
  }
  if (sem.Nested("solve")) {
    auto info = t.solver->Solve(fc_system, fc_init, fc_sol, m);
    return {info.residual, info.iter};
  }
  return {};
}

Solver GetSolver(std::string name) {
  if (name == "hypre") {
    return Solver::hypre;
  } else if (name == "zero") {
    return Solver::zero;
  } else if (name == "jacobi") {
    return Solver::jacobi;
  } else if (name == "conjugate") {
    return Solver::conjugate;
  }
  fassert(false, "Unknown solver=" + name);
}

SolverInfo Solve(
    M& m, const FieldCell<Expr>& fc_system, FieldCell<Scal>& fc_sol,
    Solver solver, int maxiter, Scal tol) {
  switch (solver) {
    case Solver::hypre:
      return SolveHypre(fc_system, &fc_sol, fc_sol, m, maxiter, tol);
    case Solver::zero:
      fc_sol.Reinit(m, 0);
      return SolverInfo{0, 0};
    case Solver::jacobi:
      return SolveJacobi(m, fc_system, &fc_sol, fc_sol, maxiter, tol);
    case Solver::conjugate:
      return SolveConjugate(m, fc_system, &fc_sol, fc_sol, maxiter, tol);
  }
  fassert(false);
}

void Run(M& m, Vars& var) {
  auto sem = m.GetSem(__func__);
  struct {
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
      t.fc_sol_exact[c] = std::sin(2 * M_PI * std::pow(x[0], 1)) *
                          std::sin(2 * M_PI * std::pow(x[1], 2)) *
                          std::sin(2 * M_PI * std::pow(x[2], 3));
    }
    m.Comm(&t.fc_sol_exact);
  }
  if (sem()) {
    // resistivity
    t.ff_rho.Reinit(m);
    for (auto f : m.FacesM()) {
      t.ff_rho[f] = (f.center().dist(m.GetGlobalLength() * 0.5) < 0.2 ? 10 : 1);
    }

    // system, only coefficients, zero constant term
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

    // constant term from exact solution
    for (auto c : m.Cells()) {
      t.fc_system[c].back() = -UEB::Eval(t.fc_system[c], c, t.fc_sol_exact, m);
    }

    // initial guess
    t.fc_sol.Reinit(m, 0);
    for (auto c : m.SuCellsM()) {
      t.fc_sol[c] =
          t.fc_sol_exact[c] * (t.fc_sol_exact[c] * 0.1 + 1);
    }

    t.solver = GetSolver(var.String["solver"]);
    m.flags.linreport = var.Int["VERBOSE"];
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
    if (var.Int("dump", 0)) {
      m.Dump(&t.fc_sol, "sol");
      m.Dump(&t.fc_sol_exact, "exact");
      m.Dump(&t.fc_diff, "diff");
    }
    if (var.Int["VERBOSE"] && m.IsRoot()) {
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
  ArgumentParser parser("Test for linear solvers.");
  parser.AddSwitch("--verbose").Help("Print solver info.");
  parser.AddVariable<std::string>("--solver", "hypre")
      .Help(
          "Linear solver to use."
          " Options are: hypre, zero, conjugate");
  parser.AddVariable<double>("--tol", 1e-3).Help("Convergence tolerance");
  parser.AddVariable<int>("--maxiter", 100).Help("Maximum iterations");
  parser.AddVariable<int>("--px", 2).Help("MPI ranks in x-direction");
  parser.AddVariable<int>("--bx", 1).Help("Blocks per rank in x-direction");
  parser.AddSwitch("--dump").Help(
      "Dump solution, exact solution, and difference");
  auto args = parser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  std::string conf = R"EOF(
set int bx 1
set int by 2
set int bz 2

set int bsx 16
set int bsy 16
set int bsz 16

set int px 2
set int py 1
set int pz 1
)EOF";

  conf += "\nset string solver " + args.String["solver"];
  conf += "\nset double hypre_symm_tol " + args.Double.GetStr("tol");
  conf += "\nset int hypre_symm_maxiter " + args.Int.GetStr("maxiter");
  conf += "\nset double tol " + args.Double.GetStr("tol");
  conf += "\nset int maxiter " + args.Int.GetStr("maxiter");
  conf += "\nset int dump " + args.Int.GetStr("dump");
  conf += "\nset int VERBOSE " + args.Int.GetStr("verbose");

  const auto px = args.Int["px"];
  const auto bx = args.Int["bx"];
  conf += "\nset int px " + std::to_string(px);
  conf += "\nset int bx " + std::to_string(bx);
  conf += "\nset int by " + std::to_string(px * bx);
  conf += "\nset int bz " + std::to_string(px * bx);

  return RunMpiBasic<M>(argc, argv, Run, conf);
}
