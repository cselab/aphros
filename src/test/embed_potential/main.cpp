// Created by Petr Karnakov on 22.01.2020
// Copyright 2020 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <utility>

#include "distr/distrbasic.h"
#include "dump/raw.h"
#include "func/init_bc.h"
#include "linear/linear.h"
#include "parse/argparse.h"
#include "solver/approx2.h"
#include "solver/approx_eb.h"
#include "solver/embed.h"
#include "solver/fluid.h"
#include "solver/reconst.h"
#include "util/distr.h"
#include "util/linear.h"

template <class M>
void CalcPotential(
    const MapEmbed<BCond<typename M::Scal>>& mebc, M& m, const Embed<M>& eb,
    FieldCell<typename M::Scal>& fcp, FieldEmbed<typename M::Scal>& feg,
    FieldCell<typename M::Scal>& fc_residual,
    std::shared_ptr<linear::Solver<M>> linsolver, std::string system_out) {
  using ExprFace = typename M::ExprFace;
  using Expr = typename M::Expr;

  auto sem = m.GetSem();
  struct {
    FieldEmbed<ExprFace> feg;
    FieldCell<Expr> fc_system;
  } * ctx(sem);
  auto& t = *ctx;

  if (sem("init")) {
    t.feg = UEmbed<M>::GradientImplicit(fcp, mebc, eb);
    t.fc_system.Reinit(m, Expr::GetUnit(0));
    for (auto c : eb.Cells()) {
      Expr sum(0);
      eb.LoopNci(c, [&](auto q) {
        const auto cf = eb.GetFace(c, q);
        eb.AppendExpr(
            sum, t.feg[cf] * eb.GetArea(cf) * eb.GetOutwardFactor(c, q), q);
      });
      t.fc_system[c] = sum;
    }
  }
  if (sem.Nested("dump_system")) {
    if (system_out.length()) {
      // dump::Raw<M>::Write(t.fc_system, system_out, m);
    }
  }
  if (sem.Nested("solve")) {
    linsolver->Solve(t.fc_system, &fcp, fcp, m);
  }
  if (sem("flux")) {
    feg.Reinit(m, 0);
    eb.LoopFaces([&](auto cf) { //
      feg[cf] = UEmbed<M>::Eval(t.feg[cf], cf, fcp, eb);
    });
    fc_residual.Reinit(m, 0);
    for (auto c : m.Cells()) {
      fc_residual[c] = UEmbed<M>::Eval(t.fc_system[c], c, fcp, eb);
    }
  }
  if (sem()) {
  }
}

template <class M>
void Run(M& m, Vars& var) {
  auto sem = m.GetSem(__func__);
  using Scal = typename M::Scal;
  using EB = Embed<M>;
  struct {
    std::unique_ptr<EB> eb;
    FieldCell<Scal> fcp;
    FieldEmbed<Scal> feg;
    FieldCell<Scal> fcdiv;
    FieldNode<Scal> fnl;
    std::shared_ptr<linear::Solver<M>> linsolver;
    typename UInitEmbedBc<M>::PlainBc bc;
    std::array<FieldCell<Scal>, 7> comps;
    FieldCell<Scal> fc_residual;
    std::string system_out;
    size_t nsteps;
    FieldCell<Scal> fc_gradexp;
  } * ctx(sem);
  auto& t = *ctx;

  if (sem("ctor")) {
    m.flags.linreport = var.Int["linreport"];
    m.flags.check_symmetry = var.Int["check_symmetry"];
    t.eb.reset(new EB(m));
    t.linsolver = ULinear<M>::MakeLinearSolver(var, "symm", m);
    t.system_out = var.String["system_out"];
    t.nsteps = var.Int["nsteps"];
    t.fcp.Reinit(m, 0);
  }
  if (sem.Nested("levelset")) {
    UEmbed<M>::InitLevelSet(t.fnl, m, var, m.IsRoot());
  }
  if (sem.Nested("init")) {
    t.eb->Init(t.fnl);
  }
  if (sem("bc")) {
    t.bc = UInitEmbedBc<M>::GetPlainBc(var.String["bc_path"], *t.eb, {});
    if (m.IsRoot()) {
      std::ofstream fdesc("bc_groups.dat");
      for (size_t i = 0; i < t.bc.vdesc.size(); ++i) {
        fdesc << i << " " << t.bc.vdesc[i] << std::endl;
      }
    }
  }
  if (sem.Nested("dump-eb")) {
    t.eb->DumpPoly();
  }
  if (sem.Nested("dump-bc")) {
    UInitEmbedBc<M>::DumpBcPoly("bc.vtk", t.bc.me_group, *t.eb, m);
  }
  for (size_t i = 0; i < t.nsteps; ++i) {
    if (sem.Nested("potential")) {
      CalcPotential(
          t.bc.mebc, m, *t.eb, t.fcp, t.feg, t.fc_residual, t.linsolver,
          t.system_out);
    }
    if (sem("dump-solution")) {
      auto& eb = *t.eb;
      t.fcdiv = Approx2<EB>::GetRegularDivergence(t.feg, *t.eb);

      if (var.Int["dump_fluxes"]) {
        for (size_t q = 0; q < 7; ++q) {
          t.comps[q].Reinit(m, 0);
        }
        for (auto c : m.Cells()) {
          for (auto q : m.Nci(c)) {
            t.comps[q.raw()][c] = t.feg[m.GetFace(c, q)];
          }
          t.comps[6][c] = t.feg[c];
        }
        m.Dump(&t.comps[0], "feg.xm");
        m.Dump(&t.comps[1], "feg.xp");
        m.Dump(&t.comps[2], "feg.ym");
        m.Dump(&t.comps[3], "feg.yp");
        m.Dump(&t.comps[4], "feg.zm");
        m.Dump(&t.comps[5], "feg.zp");
        m.Dump(&t.comps[6], "feg.c");
      }

      t.fc_gradexp.Reinit(m, 0);
      auto feg = UEmbed<M>::Gradient(t.fcp, t.bc.mebc, eb);
      for (auto c : eb.CFaces()) {
        t.fc_gradexp[c] = feg[c];
      }
      m.Dump(&t.fc_gradexp, "gradexp");

      m.Dump(&t.fc_residual, "residual");

      m.Dump(&t.fcdiv, "div");
      m.Dump(&t.fcp, "p");
    }
  }
  if (sem()) {
  }
}

template <size_t dim>
int RunGeneric(MpiWrapper& mpi, const Vars& args, std::string conf) {
  using M = MeshStructured<double, dim>;
  using MIdx = typename M::MIdx;
  MIdx meshsize(1);
  MIdx blocksize(1);
  meshsize[0] = args.Int["nx"];
  meshsize[1] = args.Int["nx"];
  blocksize[0] = args.Int["bs"];
  blocksize[1] = args.Int["bs"];
  Subdomains<MIdx> sub(meshsize, blocksize, mpi.GetCommSize());
  conf += sub.GetConfig();

  return RunMpiBasicString<M>(mpi, Run<M>, conf);
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
  ArgumentParser parser("Solver for the Poisson equation", mpi.IsRoot());
  parser.AddVariable<int>("--nx", 128).Help("Mesh size in x,y");
  parser.AddVariable<int>("--dim", 3).Options({2, 3}).Help("Space dimension");
  parser.AddVariable<int>("--nsteps", 1)
      .Help(
          "Number of steps to dump, repeated with the current solution as "
          "initial guess");
  parser.AddVariable<int>("--bs", 16)
      .Help("Block size in x,y")
      .Options({8, 16, 32});
  parser.AddSwitch("--dump_fluxes").Help("Dump fluxes from solution");
  parser.AddVariable<std::string>("--system_out", "")
      .Help(
          "Output path for linear system with one field 'data' of "
          "size (nz,ny,nx,8)");
  parser.AddVariable<std::string>("--extra", "")
      .Help("Extra configuration (commands 'set ... ')");
  auto args = parser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  std::string conf;
  conf += R"EOF(
set string backend native
set string dumpformat raw
set string eb_init list
set string eb_list_path body.dat
set string bc_path bc.dat
set int eb_init_inverse 1

set int linreport 0
set int check_symmetry 0

set int dim 2

set string linsolver_symm conjugate

set int hypre_periodic_x 0
set int hypre_periodic_y 0
set int hypre_periodic_z 1
set double hypre_symm_tol 1e-8
set int hypre_symm_maxiter 100
)EOF";

  conf += "\n" + args.String["extra"];
  conf += "\nset int nsteps " + args.Int.GetStr("nsteps");
  conf += "\nset int dump_fluxes " + args.Int.GetStr("dump_fluxes");
  conf += "\nset string system_out " + args.String.GetStr("system_out");
  conf += "\n";

  const auto dim = args.Int["dim"];
  switch (dim) {
#if USEFLAG(DIM2)
    case 2:
      return RunGeneric<2>(mpi, args, conf);
#endif
#if USEFLAG(DIM3)
    case 3:
      return RunGeneric<3>(mpi, args, conf);
#endif
    default:
      fassert(false, "Unknown dim=" + std::to_string(dim));
  }
}
