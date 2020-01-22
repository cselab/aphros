// Created by Petr Karnakov on 22.01.2020
// Copyright 2020 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <utility>

#include "distr/distrbasic.h"
#include "linear/linear.h"
#include "solver/convdiffv_eb.h"
#include "solver/embed.h"
#include "solver/fluid.h"
#include "solver/reconst.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using EB = Embed<M>;
using Type = typename EB::Type;

template <class M>
void CalcPotential(
    const MapCondFace& mfc, M& m, const Embed<M>& eb, FieldCell<Scal>& fcp,
    FieldFace<Scal>& ffv) {
  using Scal = typename M::Scal;
  using ExprFace = GVect<Scal, 3>;
  using Expr = GVect<Scal, M::dim * 2 + 2>;

  auto sem = m.GetSem();
  struct {
    FieldFace<ExprFace> ffe; // expression for volume flux in terms of p [i]
    FieldCell<Expr> fce; // linear system for potential [i]
  } * ctx(sem);
  auto& ffe = ctx->ffe;
  auto& fce = ctx->fce;

  if (sem("init")) {
    ffe.Reinit(m);
    for (auto f : eb.Faces()) {
      ExprFace e(0);
      Scal a = -eb.GetArea(f);
      e[0] = -a;
      e[1] = a;
      ffe[f] = e;
    }
    // overwrite boundary conditions
    // TODO
    // ...

    // initialize as diagonal system
    fce.Reinit(m, Expr::GetUnit(0));
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
      fce[c] = e;
    }
  }
  if (sem("solve")) {
    std::vector<Scal>* lsa;
    std::vector<Scal>* lsb;
    std::vector<Scal>* lsx;
    m.GetSolveTmp(lsa, lsb, lsx);
    lsx->resize(m.GetInBlockCells().size());
    size_t i = 0;
    for (auto c : m.Cells()) {
      (void)c;
      (*lsx)[i++] = 0;
    }
    auto l = ConvertLsCompact(fce, *lsa, *lsb, *lsx, m);
    using T = typename M::LS::T;
    l.t = T::symm;
    m.Solve(l);
  }
  if (sem("copy")) {
    std::vector<Scal>* lsa;
    std::vector<Scal>* lsb;
    std::vector<Scal>* lsx;
    m.GetSolveTmp(lsa, lsb, lsx);

    fcp.Reinit(m);
    size_t i = 0;
    for (auto c : m.Cells()) {
      fcp[c] = (*lsx)[i++];
    }
    m.Comm(&fcp);
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

void Run(M& m, Vars& var) {
  auto sem = m.GetSem("Run");
  struct {
    std::unique_ptr<EB> eb;
    MapCondFace mfc;
    FieldCell<Scal> fcp;
    FieldFace<Scal> ffv;
    FieldCell<Scal> fcdiv;
    FieldNode<Scal> fnl;
  } * ctx(sem);
  auto& fcdiv = ctx->fcdiv;
  auto& mfc = ctx->mfc;
  auto& fcp = ctx->fcp;
  auto& ffv = ctx->ffv;

  if (sem("ctor")) {
    ctx->eb.reset(new EB(m));
    ctx->fnl = InitEmbed(m, var, m.IsRoot());
  }
  if (sem.Nested("init")) {
    ctx->eb->Init(ctx->fnl);
  }
  if (sem.Nested("dumppoly")) {
    ctx->eb->DumpPoly();
  }
  if (sem.Nested("flux-proj")) {
    CalcPotential(mfc, m, *ctx->eb, fcp, ffv);
  }
  if (sem("dump")) {
    m.Dump(&fcdiv, "div");
    m.Dump(&ctx->fcp, "p");
  }
  if (sem()) {
  }
}

int main(int argc, const char** argv) {
  std::string conf = R"EOF(
set int bx 1
set int by 1
set int bz 1

set int bsx 16
set int bsy 16
set int bsz 1

set string eb_init sphere
set vect eb_sphere_c 0.5 0.5 0.5
set vect eb_sphere_r 0.249 0.249 0.249
set double eb_sphere_angle 0
set int eb_init_inverse 0

set int hypre_periodic_z 1

)EOF";
  return RunMpiBasic<M>(argc, argv, Run, conf);
}
