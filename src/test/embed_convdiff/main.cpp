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
#include "solver/approx2.h"
#include "solver/approx_eb.h"
#include "solver/convdiffe.h"
#include "solver/convdiffvg.h"
#include "solver/embed.h"
#include "solver/fluid.h"
#include "solver/reconst.h"

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using EB = Embed<M>;
using Type = typename EB::Type;
using CD = ConvDiffVectGeneric<EB, ConvDiffScalExp<EB>>;

void Run(M& m, Vars& var) {
  auto sem = m.GetSem("Run");
  struct {
    std::unique_ptr<EB> eb;
    std::unique_ptr<CD> cd;
    MapEmbed<BCond<Vect>> mebc;
    FieldCell<Scal> fcr;
    FieldCell<Scal> fcp;
    FieldEmbed<Scal> fed;
    FieldCell<Vect> fcs;
    FieldEmbed<Scal> fev;
    FieldCell<Scal> fcdiv;
    FieldNode<Scal> fnl;
    size_t frame = 0;
  } * ctx(sem);
  auto& cd = ctx->cd;
  auto& fev = ctx->fev;
  auto& fcdiv = ctx->fcdiv;
  auto& mebc = ctx->mebc;
  auto& frame = ctx->frame;

  const Vect vel(1., 1., 0.);

  if (sem("ctor")) {
    ctx->eb.reset(new EB(m));
  }
  if (sem.Nested("levelset")) {
    UEmbed<M>::InitLevelSet(ctx->fnl, m, var, m.IsRoot());
  }
  if (sem.Nested("init")) {
    ctx->eb->Init(ctx->fnl);
  }
  if (sem("init")) {
    auto& eb = *ctx->eb;
    ctx->fcr.Reinit(m, 1);
    ctx->fed.Reinit(m, 1);
    ctx->fcs.Reinit(m, Vect(0));
    ctx->fcp.Reinit(m, 0);

    fev.Reinit(m, 0);
    for (auto f : m.Faces()) {
      fev[f] = vel.dot(m.GetSurface(f));
    }

    FieldCell<Vect> fcvel(m, vel);
    const Scal dt0 = 0.5 * sqr(m.GetCellSize()[0]);
    const Scal dt0a = 1 / vel.abs().max() * m.GetCellSize()[0];
    const Scal dt = dt0 * 0.5;
    // const Scal dt = dt0a * 0.05;
    std::cout << "dt0=" << dt0 << " "
              << "dt0a=" << dt0a << " "
              << "dt=" << dt << " "
              << "dt/dt0=" << dt / dt0 << "dt/dt0a=" << dt / dt0a << " "
              << std::endl;
    typename CD::Par par;
    typename CD::Args args{fcvel, mebc, &ctx->fcr, &ctx->fed, &ctx->fcs,
                           &fev,  0,    dt,        nullptr,   par};
    cd.reset(new CD(m, eb, args));
  }
  if (sem.Nested("dumppoly")) {
    ctx->eb->DumpPoly();
  }
  const size_t maxt = 100;
  const size_t nfr = 100;
  for (size_t t = 0; t < maxt; ++t) {
    if (sem("flux-init")) {
      auto& eb = *ctx->eb;
      auto fevel = UEmbed<M>::Interpolate(
          cd->GetVelocity(), MapEmbed<BCond<Vect>>(), eb);
      for (auto f : eb.Faces()) {
        fev[f] = fevel[f].dot(eb.GetNormal(f) * eb.GetArea(f));
      }
    }
    if (sem("div")) {
      auto& eb = *ctx->eb;
      fcdiv = Approx2<EB>::GetRegularDivergence(fev, eb);
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
set double eb_box_angle 0

set vect eb_box_r 0.3 0.3 0.3

#set string eb_init sphere
set vect eb_sphere_c 0.5 0.5 0.5
set vect eb_sphere_r 0.249 0.249 0.249
set double eb_sphere_angle 0

set double cflvis 0.5

set int eb_init_inverse 0
)EOF";

  MpiWrapper mpi(&argc, &argv);
  return RunMpiBasicString<M>(mpi, Run, conf);
}
