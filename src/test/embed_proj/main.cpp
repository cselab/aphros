// Created by Petr Karnakov on 18.01.2020
// Copyright 2020 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <streambuf>
#include <string>
#include <utility>

#include "distr/distrbasic.h"
#include "linear/linear.h"
#include "solver/approx_eb.h"
#include "solver/convdiffv_eb.h"
#include "solver/embed.h"
#include "solver/fluid.h"
#include "solver/proj.h"
#include "solver/reconst.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using EB = Embed<M>;
using UEB = UEmbed<M>;
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
    std::unique_ptr<Proj<EB>> fs;
    MapCondFaceFluid mfc;
    FieldCell<Scal> fcdiv;
    FieldCell<Scal> fcr;
    FieldCell<Scal> fcd;
    FieldCell<Vect> fcf;
    FieldFace<Scal> ffbp;
    FieldCell<Scal> fcsv;
    FieldCell<Scal> fcsm;
    FieldNode<Scal> fnl;
    size_t frame = 0;
  } * ctx(sem);
  auto& fs = ctx->fs;
  auto& fcdiv = ctx->fcdiv;
  auto& mfc = ctx->mfc;
  auto& frame = ctx->frame;

  const Vect vel(var.Vect["vel"]);

  if (sem("ctor")) {
    ctx->eb.reset(new EB(m));
    ctx->fnl = UEB::InitEmbed(m, var, m.IsRoot());
  }
  if (sem.Nested("init")) {
    ctx->eb->Init(ctx->fnl);
  }
  if (sem("init")) {
    auto& eb = *ctx->eb;
    ctx->fcr.Reinit(m, 1);
    ctx->fcd.Reinit(m, 1);
    ctx->fcf.Reinit(m, Vect(0));
    ctx->fcsv.Reinit(m, 0);
    ctx->fcsm.Reinit(m, 0);
    ctx->ffbp.Reinit(m, 0);

    const Vect force(var.Vect["force"]);
    for (auto f : m.AllFaces()) {
      ctx->ffbp[f] = force.dot(m.GetNormal(f));
    }

    FieldCell<Vect> fcvel(m, vel);
    const Scal dt0 = 0.5 * sqr(m.GetCellSize()[0]);
    const Scal dt0a = 1 / vel.abs().max() * m.GetCellSize()[0];
    const Scal cflvis = var.Double["cflvis"];
    const Scal dt = dt0 * cflvis;
    // const Scal dt = dt0a * 0.05;
    std::cout << "dt0=" << dt0 << " "
              << "dt0a=" << dt0a << " "
              << "dt=" << dt << " "
              << "dt/dt0=" << dt / dt0 << "dt/dt0a=" << dt / dt0a << " "
              << std::endl;
    typename Proj<EB>::Par par;
    fs.reset(new Proj<EB>(
        m, eb, fcvel, mfc, MapCell<std::shared_ptr<CondCellFluid>>(), &ctx->fcr,
        &ctx->fcd, &ctx->fcf, &ctx->ffbp, &ctx->fcsv, &ctx->fcsm, 0, dt, par));
  }
  if (sem.Nested("dumppoly")) {
    ctx->eb->DumpPoly();
  }
  const size_t maxt = 100;
  const size_t nfr = 100;
  for (size_t t = 0; t < maxt; ++t) {
    if (sem("div")) {
      auto& eb = *ctx->eb;
      fcdiv = GetDivergence(fs->GetVolumeFlux(), m, eb);
    }
    if (sem.Nested("convdiff-start")) {
      fs->StartStep();
    }
    if (sem.Nested("convdiff-iter")) {
      fs->MakeIteration();
    }
    if (sem.Nested("convdiff-finish")) {
      fs->FinishStep();
    }
    if (t % std::max<size_t>(1, maxt / nfr) != 0) {
      continue;
    }
    if (sem("dump")) {
      m.Dump(&fs->GetVelocity(), 0, "vx");
      m.Dump(&fs->GetVelocity(), 1, "vy");
      m.Dump(&fs->GetVelocity(), 2, "vz");
      m.Dump(&fs->GetPressure(), "p");
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
set vect eb_box_r 2 0.251 2
#set vect eb_box_r 10 0.249 10
set double eb_box_angle 0.05

set string eb_init sphere
set vect eb_sphere_c 0.5 0.5 0.5
set vect eb_sphere_r 0.3 0.2 0.2
set double eb_sphere_angle 0.2

set vect vel 0 0 0
set vect force 1 0 0
set double cflvis 0.5

set int eb_init_inverse 0
)EOF";

  std::ifstream fin("add.conf");
  std::string add;
  add.assign(
      std::istreambuf_iterator<char>(fin), std::istreambuf_iterator<char>());
  return RunMpiBasic<M>(argc, argv, Run, conf + add);
}
