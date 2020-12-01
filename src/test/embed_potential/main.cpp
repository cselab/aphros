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
#include "func/init_bc.h"
#include "linear/linear.h"
#include "solver/approx2.h"
#include "solver/approx_eb.h"
#include "solver/convdiffv_eb.h"
#include "solver/embed.h"
#include "solver/fluid.h"
#include "solver/reconst.h"
#include "util/linear.h"

template <class M>
void CalcPotential(
    const MapEmbed<BCond<typename M::Scal>>& mebc, M& m, const Embed<M>& eb,
    FieldCell<typename M::Scal>& fcp, FieldEmbed<typename M::Scal>& ffv,
    std::shared_ptr<linear::Solver<M>> linsolver) {
  using Scal = typename M::Scal;
  using ExprFace = generic::Vect<Scal, 3>;
  using Expr = generic::Vect<Scal, M::dim * 2 + 2>;

  auto sem = m.GetSem();
  struct {
    FieldFace<ExprFace> ffe; // expression for volume flux in terms of p [i]
    FieldCell<Expr> fc_system; // linear system for potential [i]
  } * ctx(sem);
  auto& ffe = ctx->ffe;
  auto& fc_system = ctx->fc_system;

  if (sem("init")) {
    ffe.Reinit(m);
    for (auto f : eb.Faces()) {
      ExprFace e(0);
      const IdxCell cm = m.GetCell(f, 0);
      const IdxCell cp = m.GetCell(f, 1);
      Scal dn =
          (eb.GetCellCenter(cp) - eb.GetCellCenter(cm)).dot(eb.GetNormal(f));
      dn = (dn > 0 ? 1 : -1) * std::max(std::abs(dn), m.GetCellSize()[0] * 0.5);
      const Scal a = eb.GetArea(f) / dn;
      e[0] = -a;
      e[1] = a;
      ffe[f] = e;
    }
    // overwrite boundary conditions
    for (auto& it : mebc.GetMapFace()) {
      const IdxFace f = it.first;
      const auto& bc = it.second;
      if (bc.type == BCondType::neumann) {
        ExprFace e(0);
        e[2] = bc.val * eb.GetArea(f);
        ffe[f] = e;
      } else {
        throw std::runtime_error("unknown bc");
      }
    }

    // initialize as diagonal system
    fc_system.Reinit(m, Expr::GetUnit(0));
    // overwrite with div=0 equation in non-excluded cells
    for (auto c : eb.Cells()) {
      Expr e(0);
      for (auto q : eb.Nci(c)) {
        const IdxFace f = m.GetFace(c, q);
        const ExprFace v = ffe[f] * m.GetOutwardFactor(c, q);
        e[0] += v[1 - q % 2];
        e[1 + q] += v[q % 2];
        e[Expr::dim - 1] += v[2];
      }
      fc_system[c] = e;
    }
  }
  if (sem.Nested("solve")) {
    linsolver->Solve(fc_system, nullptr, fcp, m);
  }
  if (sem("flux")) {
    ffv.Reinit(m);
    for (auto f : eb.Faces()) {
      const IdxCell cm = m.GetCell(f, 0);
      const IdxCell cp = m.GetCell(f, 1);
      const auto& e = ffe[f];
      ffv[f] = e[0] * fcp[cm] + e[1] * fcp[cp] + e[2];
    }
  }
}

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using EB = Embed<M>;
using UEB = UEmbed<M>;
using Type = typename EB::Type;

void Run(M& m, Vars& var) {
  auto sem = m.GetSem("Run");
  struct {
    std::unique_ptr<EB> eb;
    MapEmbed<BCond<Scal>> mebc;
    FieldCell<Scal> fcp;
    FieldEmbed<Scal> ffv;
    FieldCell<Scal> fcdiv;
    FieldNode<Scal> fnl;
    std::shared_ptr<linear::Solver<M>> linsolver;
  } * ctx(sem);
  auto& t = *ctx;

  if (sem("eb_ctor")) {
    t.eb.reset(new EB(m));
    t.linsolver = ULinear<M>::MakeLinearSolver(var, "symm", m);
  }
  if (sem.Nested("levelset")) {
    UEB::InitLevelSet(t.fnl, m, var, m.IsRoot());
  }
  if (sem.Nested("eb_init")) {
    t.eb->Init(t.fnl);
  }
  if (sem("eb_ctor")) {
    auto& eb = *t.eb;
    auto res = UInitEmbedBc<M>::GetPlainBc(var.String["bc_path"], eb, {});
    t.mebc = res.mebc;
  }
  if (sem.Nested("eb_dumppoly")) {
    t.eb->DumpPoly();
  }
  if (sem.Nested("flux-proj")) {
    CalcPotential(t.mebc, m, *t.eb, t.fcp, t.ffv, t.linsolver);
  }
  if (sem("dump")) {
    t.fcdiv = Approx2<EB>::GetRegularDivergence(t.ffv, *t.eb);
    m.Dump(&t.fcdiv, "div");
    m.Dump(&t.fcp, "p");
  }
  if (sem()) {
  }
}

int main(int argc, const char** argv) {
  const int bsx = 32;
  int nx = 128; // mesh size
  if (argc > 1) {
    nx = atoi(argv[1]);
  }
  assert(nx % bsx == 0);
  const int bx = nx / bsx;

  std::string conf = R"EOF(
set int bz 1
set int bsz 1

set string eb_init list
set string eb_list_path body.dat
set int eb_init_inverse 1

set string dumpformat plain

set int dim 2

set int hypre_periodic_z 1
set double hypre_symm_tol 1e-8
set int hypre_symm_maxiter 10000
)EOF";

  conf += "set int bx " + std::to_string(bx) + "\n";
  conf += "set int by " + std::to_string(bx) + "\n";
  conf += "set int bsx " + std::to_string(bsx) + "\n";
  conf += "set int bsy " + std::to_string(bsx) + "\n";

  MpiWrapper mpi(&argc, &argv);
  return RunMpiBasicString<M>(mpi, Run, conf);
}
