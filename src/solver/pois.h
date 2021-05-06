// Created by Petr Karnakov on 06.11.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <cmath>
#include <exception>
#include <memory>
#include <sstream>
#include <string>

#include "cond.h"
#include "geom/mesh.h"
#include "linear/linear.h"
#include "solver/approx_eb.h"

// Solves Poisson equation: laplace u = rhs.
// fc_rhs: rhs [i]
// mebc: boundary conditions
// centered_rhs: subtract average from the rhs to make the problem
//               with zero-Neumann conditions consistent
// Output:
// fc_sol: solution [a]
template <class M>
void SolvePoisson(
    FieldCell<typename M::Scal>& fcu, const FieldCell<typename M::Scal>& fc_rhs,
    const MapEmbed<BCond<typename M::Scal>>& mebc,
    std::shared_ptr<linear::Solver<M>> linsolver, M& m,
    bool centered_rhs=true) {
  using Scal = typename M::Scal;
  using Expr = typename M::Expr;
  using ExprFace = typename M::ExprFace;
  auto sem = m.GetSem("pois");
  struct {
    Scal avg_rhs = 0;
    Scal sum_vol = 0;
    FieldCell<Expr> fcl;
  } * ctx(sem);
  auto& t = *ctx;
  if (sem("reduce")) {
    for (auto c : m.Cells()) {
      auto v = m.GetVolume(c);
      t.avg_rhs += fc_rhs[c] * v;
      t.sum_vol += v;
    }
    m.Reduce(&t.avg_rhs, "sum");
    m.Reduce(&t.sum_vol, "sum");
  }
  if (sem("assemble")) {
    t.avg_rhs /= t.sum_vol;
    const auto ffg = UEmbed<M>::GradientImplicit(mebc, m);
    t.fcl.Reinit(m, Expr::GetUnit(0));
    for (auto c : m.Cells()) {
      Expr sum(0);
      m.LoopNci(c, [&](auto q) {
        const auto cf = m.GetFace(c, q);
        const ExprFace flux = ffg[cf] * m.GetArea(cf);
        m.AppendExpr(sum, flux * m.GetOutwardFactor(c, q), q);
      });
      const auto rhs = fc_rhs[c] - (centered_rhs ? t.avg_rhs : 0);
      sum.back() = -rhs * m.GetVolume(c);
      t.fcl[c] = sum;
    }
  }
  if (sem.Nested("solve")) {
    linsolver->Solve(t.fcl, nullptr, fcu, m);
  }
  if (sem()) {
  }
}
