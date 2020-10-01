// Created by Petr Karnakov on 01.10.2020
// Copyright 2020 ETH Zurich

#include <iostream>

#include "distr/distrbasic.h"
#include "solver/approx_eb.h"
#include "solver/embed.h"
#include "linear/linear.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

void Run(M& m, Vars&) {
  auto sem = m.GetSem();
  using Expr = typename M::Expr;
  using ExprFace = typename M::ExprFace;
  struct {
    FieldCell<Scal> fc_rhs;
    FieldCell<Scal> fcu;
    FieldFace<Scal> ff_rho;
    FieldCell<Expr> fc_system;
    MapEmbed<BCond<Scal>> mebc;
  } * ctx(sem);
  auto& t = *ctx;
  if (sem()) {
    // rhs
    t.fc_rhs.Reinit(m, 0);
    for (auto c : m.CellsM()) {
      t.fc_rhs[c] = c.center[0];
    }

    // resistivity
    t.ff_rho.Reinit(m);
    for (auto f : m.FacesM()) {
      t.ff_rho[f] = (Vect(0.5, 0.25, 0.125).dist(f.center()) < 0.2 ? 1000 : 1);
    }

    // initial guess
    t.fcu.Reinit(m, 0);

    // system
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
  }
  if (sem.Nested("solve")) {
    Solve(t.fc_system, &t.fcu, t.fcu, M::LS::T::symm, m, "symm");
  }
}

int main(int argc, const char** argv) {
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

  return RunMpiBasic<M>(argc, argv, Run, conf);
}
