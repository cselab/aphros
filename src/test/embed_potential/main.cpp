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
#include "linear/linear.h"
#include "solver/convdiffv_eb.h"
#include "solver/embed.h"
#include "solver/approx_eb.h"
#include "solver/fluid.h"
#include "solver/reconst.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using EB = Embed<M>;
using UEB = UEmbed<M>;
using Type = typename EB::Type;

template <class M>
FieldCell<typename M::Scal> GetDivergence(
    FieldFace<typename M::Scal>& ffv, const Embed<M>& eb) {
  auto& m = eb.GetMesh();
  using Scal = typename M::Scal;
  FieldCell<Scal> fcdiv(m, 0);
  for (auto c : eb.Cells()) {
    Scal div = 0;
    for (auto q : eb.Nci(c)) {
      div += ffv[m.GetFace(c, q)] * m.GetOutwardFactor(c, q);
    }
    fcdiv[c] = div / eb.GetVolume(c);
  }
  return fcdiv;
}

template <class M>
void CalcPotential(
    const MapCondFace& mfc, M& m, const Embed<M>& eb, FieldCell<Scal>& fcp,
    FieldFace<Scal>& ffv) {
  using Scal = typename M::Scal;
  using ExprFace = generic::Vect<Scal, 3>;
  using Expr = generic::Vect<Scal, M::dim * 2 + 2>;

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
    for (auto& it : mfc) {
      const IdxFace f = it.first;
      const auto& cb = it.second;
      if (auto* cd = cb.template Get<CondFaceGrad<Scal>>()) {
        ExprFace e(0);
        e[2] = cd->GetGrad() * eb.GetArea(f);
        ffe[f] = e;
      } else {
        throw std::runtime_error("unknown bc");
      }
    }

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
  auto& eb_ = ctx->eb;

  if (sem("eb_ctor")) {
    eb_.reset(new EB(m));
    ctx->fnl = UEB::InitEmbed(m, var, m.IsRoot());
  }
  if (sem.Nested("eb_init")) {
    eb_->Init(ctx->fnl);
  }
  if (sem("eb_ctor")) {
    auto& eb = *eb_;
    auto isclose = [](Scal a, Scal b) { return std::abs(a - b) < 1e-6; };
    for (auto f : eb.Faces()) {
      const Vect x = eb.GetFaceCenter(f);
      if (isclose(x[0], 0.)) {
        mfc[f] = UniquePtr<CondFaceGradFixed<Scal>>(1., 1);
      }
      if (isclose(x[0], 1.)) {
        mfc[f] = UniquePtr<CondFaceGradFixed<Scal>>(1., 0);
      }
      if (isclose(x[1], 0.)) {
        mfc[f] = UniquePtr<CondFaceGradFixed<Scal>>(0., 1);
      }
      if (isclose(x[1], 1.)) {
        mfc[f] = UniquePtr<CondFaceGradFixed<Scal>>(0., 0);
      }
    }
  }
  if (sem.Nested("eb_dumppoly")) {
    eb_->DumpPoly();
  }
  if (sem.Nested("flux-proj")) {
    CalcPotential(mfc, m, *eb_, fcp, ffv);
  }
  if (sem("dump")) {
    fcdiv = GetDivergence(ffv, *eb_);
    m.Dump(&fcdiv, "div");
    m.Dump(&ctx->fcp, "p");
  }
  if (sem()) {
    if (m.IsRoot()) {
      std::cout << "iter=" << m.GetIter() << std::endl;
    }
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

  return RunMpiBasic<M>(argc, argv, Run, conf);
}
