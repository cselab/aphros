// Created by Petr Karnakov on 24.12.2019
// Copyright 2019 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <iostream>
#include <memory>
#include <string>

#include <distr/distrbasic.h>
#include <solver/vof.h>

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;

void Run(M& m, Vars&) {
  auto sem = m.GetSem();

  struct {
    std::unique_ptr<Vof<M>> as; // advection solver
    FieldCell<Scal> fc_src; // volume source
    FieldEmbed<Scal> fe_flux; // volume flux
    MapEmbed<BCondAdvection<Scal>> mf_cond; // face conditions
  } * ctx(sem);

  auto& as = ctx->as;
  const Scal tmax = 1.;
  const Vect vel(0.5, 0.3, 0.1);
  const Scal cfl = 0.5;

  if (sem("init")) {
    auto& fc_src = ctx->fc_src;
    auto& fe_flux = ctx->fe_flux;
    auto& mf_cond = ctx->mf_cond;

    fc_src.Reinit(m, 0);
    fe_flux.Reinit(m, 0);
    for (auto f : m.Faces()) {
      fe_flux[f] = vel.dot(m.GetSurface(f));
    }
    FieldCell<Scal> fccl(m, 0); // initial color
    FieldCell<Scal> fcu(m, 0); // initial volume fraction
    for (auto c : m.Cells()) {
      fcu[c] = (m.GetCenter(c).dist(Vect(0.5, 0.5, 0.5)) < 0.2);
    }
    const Scal dt = cfl * m.GetCellSize()[0] / vel.norm();
    typename Vof<M>::Par p;
    as.reset(
        new Vof<M>(m, m, fcu, fccl, mf_cond, &fe_flux, &fc_src, 0., dt, p));
  }
  sem.LoopBegin();
  if (sem.Nested("start")) {
    as->StartStep();
  }
  if (sem.Nested("iter")) {
    as->MakeIteration();
  }
  if (sem.Nested("finish")) {
    as->FinishStep();
  }
  if (sem("checkloop")) {
    if (as->GetTime() >= tmax) {
      sem.LoopBreak();
    }
  }
  if (sem("dump")) {
    m.Dump(&as->GetField(), "u");
  }
  if (sem()) {
  }
  sem.LoopEnd();
}

int main(int argc, const char** argv) {
  std::string conf = R"EOF(
set int bx 1
set int by 1
set int bz 1

set int bsx 16
set int bsy 16
set int bsz 16

set int px 1
set int py 1
set int pz 1
)EOF";

  MpiWrapper mpi(&argc, &argv);
  return RunMpiBasicString<M>(mpi, Run, conf);
}
