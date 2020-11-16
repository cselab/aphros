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
#include "util/distr.h"
#include "util/timer.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;
using Expr = typename M::Expr;
using ExprFace = typename M::ExprFace;
using UEB = UEmbed<M>;

void Run(M& m, Vars& var) {
  auto sem = m.GetSem(__func__);
  struct {
    FieldCell<Scal> fc_sol;
    FieldCell<Scal> fc_sol_exact;
    FieldCell<Scal> fc_diff;
    FieldFace<Scal> ff_rho;
    FieldCell<Expr> fc_system;
    MapEmbed<BCond<Scal>> mebc;
    std::unique_ptr<linear::Solver<M>> solver;
    std::vector<generic::Vect<Scal, 3>> norms;
    typename linear::Solver<M>::Info info;
    SingleTimer timer;
    double time_start;
    double time_stop;
    Scal mean_diff;
    Scal sumw;
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
      t.fc_sol[c] = t.fc_sol_exact[c] * (t.fc_sol_exact[c] * 0.1 + 1);
    }

    const auto name = var.String["linsolver_symm"];
    auto factory = linear::ModuleLinear<M>::GetInstance(name);
    fassert(factory, "Solver not found: " + name);
    t.solver = factory->Make(var, "symm", m);
    m.flags.linreport = var.Int["VERBOSE"];
    t.time_start = t.timer.GetSeconds();
  }
  if (sem.Nested("solve")) {
    t.info = t.solver->Solve(t.fc_system, &t.fc_sol, t.fc_sol, m);
  }
  if (sem("diff")) {
    t.time_stop = t.timer.GetSeconds();
    t.fc_diff.Reinit(m);
    t.mean_diff = 0;
    t.sumw = 0;
    for (auto c : m.Cells()) {
      t.fc_diff[c] = t.fc_sol[c] - t.fc_sol_exact[c];
      t.mean_diff += t.fc_diff[c];
      t.sumw += 1;
    }
    m.Reduce(&t.mean_diff, Reduction::sum);
    m.Reduce(&t.sumw, Reduction::sum);
  }
  if (sem("diff")) {
    t.mean_diff /= t.sumw;
    for (auto c : m.Cells()) {
      t.fc_diff[c] -= t.mean_diff;
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
      std::cout << "\ntime=" << std::fixed << t.time_stop - t.time_start;
      std::cout << std::endl;
    }
  }
  if (sem()) {
  }
}

int main(int argc, const char** argv) {
#if USEFLAG(HYPRE)
  FORCE_LINK(linear_hypre);
#endif
#if USEFLAG(AMGX)
  FORCE_LINK(linear_amgx);
#endif
  FORCE_LINK(linear_conjugate);
  FORCE_LINK(linear_jacobi);

  MpiWrapper mpi(&argc, &argv);

  ArgumentParser parser("Test for linear solvers.", mpi.IsRoot());
  parser.AddSwitch("--verbose").Help("Print solver info.");
  auto instances = []() {
    auto map = linear::ModuleLinear<M>::GetInstances();
    std::vector<std::string> res;
    for (auto p : map) {
      res.push_back(p.first);
    }
    return res;
  }();
  parser.AddVariable<std::string>("--solver", "hypre")
      .Help("Linear solver to use")
      .Options(instances);
  parser.AddVariable<double>("--tol", 1e-3).Help("Convergence tolerance");
  parser.AddVariable<int>("--maxiter", 100).Help("Maximum iterations");
  parser.AddVariable<int>("--mesh", 32).Help("Mesh size in all directions");
  parser.AddVariable<int>("--block", 16)
      .Help("Block size in all directions")
      .Options({8, 16, 32});
  parser.AddSwitch("--dump").Help(
      "Dump solution, exact solution, and difference");
  parser.AddVariable<std::string>("--extra", "")
      .Help("Extra configuration (commands 'set ... ')");
  auto args = parser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  Subdomains<MIdx> sub(
      MIdx(args.Int["mesh"]), MIdx(args.Int["block"]), mpi.GetCommSize());
  std::string conf = "";
  conf += sub.GetConfig();

  conf += "\nset string linsolver_symm " + args.String["solver"];
  conf += "\nset double hypre_symm_tol " + args.Double.GetStr("tol");
  conf += "\nset int hypre_symm_maxiter " + args.Int.GetStr("maxiter");
  conf += "\nset int hypre_print 0";
  conf += "\nset string hypre_symm_solver pcg";
  conf += "\nset double tol " + args.Double.GetStr("tol");
  conf += "\nset int maxiter " + args.Int.GetStr("maxiter");
  conf += "\nset int dump " + args.Int.GetStr("dump");
  conf += "\nset int VERBOSE " + args.Int.GetStr("verbose");

  conf += "\n" + args.String["extra"];

  return RunMpiBasicString<M>(mpi, Run, conf);
}
