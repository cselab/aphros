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
#include "solver/convdiffv_eb.h"
#include "solver/embed.h"
#include "solver/reconst.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using EB = Embed<M>;
using Type = typename EB::Type;
using CD = ConvDiffVectEmbed<M>;

void Run(M& m, Vars& var) {
  auto sem = m.GetSem("Run");
  struct {
    std::unique_ptr<EB> eb;
    std::unique_ptr<CD> cd;
    FieldCell<Vect> fcvel;
    MapCondFace mfc;
    FieldCell<Scal> fcr;
    FieldEmbed<Scal> fed;
    FieldCell<Vect> fcs;
    FieldFace<Scal> ffv;
    size_t frame = 0;
  } * ctx(sem);
  auto& eb = ctx->eb;
  auto& cd = ctx->cd;
  auto& fcvel = ctx->fcvel;
  auto& ffv = ctx->ffv;
  auto& mfc = ctx->mfc;
  auto& frame = ctx->frame;

  if (sem("initeb")) {
    auto fnl = InitEmbed(m, var);
    ctx->eb.reset(new EB(m, fnl));
  }
  if (sem("init")) {
    ctx->fcr.Reinit(m, 1);
    ctx->fed.Reinit(m, 0);
    ctx->fcs.Reinit(m, Vect(0));
    ffv.Reinit(m, 0);
    const Vect vel(1., 0.5, 0.25);
    for (auto f : m.Faces()) {
      ffv[f] = vel.dot(m.GetSurface(f));
    }
    fcvel.Reinit(m, Vect(1., 1., 1.));
    /*
    for (auto c : eb->AllCells()) {
      const Scal a = 12;
      fcvel[c] =
          std::sin(m.GetCenter(c)[0] * a) * std::sin(m.GetCenter(c)[1] * a);
    }
    */
    const size_t bc = 0;
    const Vect bcvel(1);
    const Scal dt0 = 0.5 * sqr(m.GetCellSize()[0]);
    const Scal dt0a = 1 / vel.abs().max() * m.GetCellSize()[0];
    //const Scal dt = dt0 * 0.5;
    const Scal dt = dt0a * 0.05;
    std::cout << "dt0=" << dt0 << " "
              << "dt0a=" << dt0a << " "
              << "dt=" << dt << " "
              << "dt/dt0=" << dt / dt0 << "dt/dt0a=" << dt / dt0a << " "
              << std::endl;
    typename CD::Par par;
    cd.reset(new CD(
        m, *eb, fcvel, mfc, bc, bcvel, &ctx->fcr, &ctx->fed, &ctx->fcs,
        &ctx->ffv, 0, dt, par));
  }
  if (sem.Nested("dumppoly")) {
    eb->DumpPoly();
  }
  const size_t maxt = 100;
  const size_t nfr = 100;
  for (size_t t = 0; t < maxt; ++t) {
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

set string eb_init sphere
set vect eb_sphere_c 0.5 0.5 0.5
set vect eb_sphere_r 0.249 0.249 0.249
set double eb_sphere_angle 0

set int eb_init_inverse 0
)EOF";
  return RunMpiBasic<M>(argc, argv, Run, conf);
}
