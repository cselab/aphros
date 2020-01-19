// Created by Petr Karnakov on 31.12.2019
// Copyright 2019 ETH Zurich

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
using CD = ConvDiffVectEmbed<M>;

template <class M>
FieldCell<typename M::Scal> GetDivergence(
    FieldEmbed<typename M::Scal>& fev, const M& m, const Embed<M>& eb) {
  using Scal = typename M::Scal;
  FieldCell<Scal> fcdiv(m, 0);
  for (auto c : eb.Cells()) {
    Scal div = fev[c];
    for (auto q : m.Nci(c)) {
      div += fev[m.GetFace(c, q)] * m.GetOutwardFactor(c, q);
    }
    fcdiv[c] = div / eb.GetVolume(c);
  }
  return fcdiv;
}

template <class M>
FieldCell<typename M::Scal> GetDivergence(
    FieldFace<typename M::Scal>& ffv, const M& m, const Embed<M>& eb) {
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
FieldCell<typename M::Scal> ProjectVolumeFlux(
    FieldFace<typename M::Scal>& ffv, const MapCondFaceFluid& mfc, M& m,
    const Embed<M>& eb) {
  using Scal = typename M::Scal;
  using ExprFace = GVect<Scal, 3>;
  using Expr = GVect<Scal, M::dim * 2 + 2>;

  auto sem = m.GetSem();
  struct {
    FieldFace<ExprFace> ffe; // expression for corrected volume flux [i]
    FieldCell<Expr> fce; // linear system for pressure [i]
    FieldFace<bool> ffbd; // true for faces with boundary conditions
    FieldCell<Scal> fcp; // pressure (up to a constant)
  } * ctx(sem);
  auto& ffe = ctx->ffe;
  auto& fce = ctx->fce;
  auto& ffbd = ctx->ffbd;
  auto& fcp = ctx->fcp;

  if (sem("init")) {
    ffbd.Reinit(m, false);
    for (auto p : mfc) {
      ffbd[p.GetIdx()] = true;
    }

    ffe.Reinit(m);
    for (auto f : eb.Faces()) {
      auto& e = ffe[f];
      if (!ffbd[f]) { // inner
        Scal a = -eb.GetArea(f);
        e[0] = -a;
        e[1] = a;
      } else { // boundary
        e[0] = 0;
        e[1] = 0;
      }
      e[2] = ffv[f];
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
  if (sem("apply")) {
    for (auto f : eb.Faces()) {
      const IdxCell cm = m.GetCell(f, 0);
      const IdxCell cp = m.GetCell(f, 1);
      const auto& e = ffe[f];
      ffv[f] = e[0] * fcp[cm] + e[1] * fcp[cp] + e[2];
    }
  }
  return fcp;
}

void Run(M& m, Vars& var) {
  auto sem = m.GetSem("Run");
  struct {
    std::unique_ptr<EB> eb;
    std::unique_ptr<CD> cd;
    MapCondFace mfc;
    FieldCell<Scal> fcr;
    FieldCell<Scal> fcp;
    FieldEmbed<Scal> fed;
    FieldCell<Vect> fcs;
    FieldFace<Scal> ffv;
    FieldCell<Scal> fcdiv;
    size_t frame = 0;
  } * ctx(sem);
  auto& cd = ctx->cd;
  auto& ffv = ctx->ffv;
  auto& fcdiv = ctx->fcdiv;
  auto& mfc = ctx->mfc;
  auto& frame = ctx->frame;

  const Vect vel(1., 1., 0.);
  const size_t bc = 0;
  const Vect bcvel(0);

  if (sem("initeb")) {
    auto fnl = InitEmbed(m, var);
    ctx->eb.reset(new EB(m, fnl));
  }
  if (sem("init")) {
    auto& eb = *ctx->eb;
    ctx->fcr.Reinit(m, 1);
    ctx->fed.Reinit(m, 1);
    ctx->fcs.Reinit(m, Vect(0));
    ctx->fcp.Reinit(m, 0);

    ffv.Reinit(m, 0);
    for (auto f : m.Faces()) {
      ffv[f] = vel.dot(m.GetSurface(f));
    }

    FieldCell<Vect> fcvel(m, vel);
    const Scal dt0 = 0.5 * sqr(m.GetCellSize()[0]);
    const Scal dt0a = 1 / vel.abs().max() * m.GetCellSize()[0];
    const Scal dt = dt0 * 0.5;
    //const Scal dt = dt0a * 0.05;
    std::cout << "dt0=" << dt0 << " "
              << "dt0a=" << dt0a << " "
              << "dt=" << dt << " "
              << "dt/dt0=" << dt / dt0 << "dt/dt0a=" << dt / dt0a << " "
              << std::endl;
    typename CD::Par par;
    cd.reset(new CD(
        m, eb, fcvel, mfc, bc, bcvel, &ctx->fcr, &ctx->fed, &ctx->fcs, &ffv, 0,
        dt, par));
  }
  if (sem.Nested("dumppoly")) {
    ctx->eb->DumpPoly();
  }
  const size_t maxt = 100;
  const size_t nfr = 100;
  for (size_t t = 0; t < maxt; ++t) {
    if (sem("flux-init")) {
      auto& eb = *ctx->eb;
      auto fevel = eb.Interpolate(cd->GetVelocity(), bc, bcvel);
      for (auto f : eb.Faces()) {
        ffv[f] = fevel[f].dot(eb.GetNormal(f) * eb.GetArea(f));
      }
    }
    if (sem.Nested("flux-proj")) {
      ctx->fcp = ProjectVolumeFlux(ffv, MapCondFaceFluid(), m, *ctx->eb);
    }
    if (sem("div")) {
      auto& eb = *ctx->eb;
      fcdiv = GetDivergence(ffv, m, eb);
    }
    if (sem.Nested("convdiff-start")) {
      cd->StartStep();
    }
    if (sem.Nested("convdiff-iter")) {
      cd->MakeIteration();
    }
    if (sem.Nested("convdiff-finish")) {
      cd->FinishStep();
    }
    if (t % std::max<size_t>(1, maxt / nfr) != 0) {
      continue;
    }
    if (sem("dump")) {
      m.Dump(&cd->GetVelocity(), 0, "vx");
      m.Dump(&cd->GetVelocity(), 1, "vy");
      m.Dump(&cd->GetVelocity(), 2, "vz");
      m.Dump(&fcdiv, "div");
      m.Dump(&ctx->fcp, "p");
      ++frame;
    }
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
set int bsz 16

set string eb_init box
set vect eb_box_c 0.5 0.5 0.5
#set vect eb_box_r 0.249 0.249 0.249
set vect eb_box_r 0.251 0.251 0.251
#set vect eb_box_r 10 0.249 10
set double eb_box_angle 0

set vect eb_box_r 0.3 0.3 0.3

#set string eb_init sphere
set vect eb_sphere_c 0.5 0.5 0.5
set vect eb_sphere_r 0.249 0.249 0.249
set double eb_sphere_angle 0

set double cflvis 0.5

set int eb_init_inverse 0
)EOF";
  return RunMpiBasic<M>(argc, argv, Run, conf);
}
