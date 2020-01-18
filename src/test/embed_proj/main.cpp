// Created by Petr Karnakov on 18.01.2020
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
#include "solver/proj_eb.h"
#include "solver/reconst.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using EB = Embed<M>;
using Type = typename EB::Type;

template <class M>
FieldCell<typename M::Scal> GetDivergence(
    const FieldEmbed<typename M::Scal>& fev, const M& m, const Embed<M>& eb) {
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
    const FieldFace<typename M::Scal>& ffv, const M& m, const Embed<M>& eb) {
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

void Run(M& m, Vars& var) {
  auto sem = m.GetSem("Run");
  struct {
    std::unique_ptr<EB> eb;
    std::unique_ptr<ProjEmbed<M>> cd;
    MapCondFaceFluid mfc;
    FieldCell<Scal> fcdiv;
    FieldCell<Scal> fcr;
    FieldCell<Scal> fcd;
    FieldCell<Vect> fcf;
    FieldFace<Scal> ffbp;
    FieldCell<Scal> fcsv;
    FieldCell<Scal> fcsm;
    size_t frame = 0;
  } * ctx(sem);
  auto& cd = ctx->cd;
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
    ctx->fcd.Reinit(m, 1);
    ctx->fcf.Reinit(m, Vect(0));
    ctx->fcsv.Reinit(m, 0);
    ctx->fcsm.Reinit(m, 0);
    ctx->ffbp.Reinit(m, 0);

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
    typename ProjEmbed<M>::Par par;
    cd.reset(new ProjEmbed<M>(
        m, fcvel, mfc, eb, bc, bcvel, MapCell<std::shared_ptr<CondCellFluid>>(),
        &ctx->fcr, &ctx->fcd, &ctx->fcf, &ctx->ffbp, &ctx->fcsv, &ctx->fcsm, 0,
        dt, par));
  }
  if (sem.Nested("dumppoly")) {
    ctx->eb->DumpPoly();
  }
  const size_t maxt = 100;
  const size_t nfr = 100;
  for (size_t t = 0; t < maxt; ++t) {
    if (sem("div")) {
      auto& eb = *ctx->eb;
      fcdiv = GetDivergence(cd->GetVolumeFlux(), m, eb);
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

set vect eb_box_r 0.3 0.3 0.3

#set string eb_init sphere
set vect eb_sphere_c 0.5 0.5 0.5
set vect eb_sphere_r 0.249 0.249 0.249
set double eb_sphere_angle 0

set int eb_init_inverse 0
)EOF";
  return RunMpiBasic<M>(argc, argv, Run, conf);
}
