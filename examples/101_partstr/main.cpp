// Created by Petr Karnakov on 25.12.2019
// Copyright 2019 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <iostream>
#include <memory>
#include <string>

#include <distr/distrbasic.h>
#include <solver/curv.h>

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;

void Run(M& m, Vars&) {
  using PSM = PartStrMeshM<M>;
  using Plic = generic::Plic<Vect>;
  auto sem = m.GetSem();

  struct {
    FieldCell<Scal> fcu; // volume fraction
    FieldCell<Scal> fca; // plane constant
    FieldCell<Vect> fcn; // normal
    FieldCell<bool> fci; // interface mask
    FieldCell<Scal> fck; // curvature
    MapEmbed<BCondAdvection<Scal>> mebc;
    typename PSM::Par par;
    std::unique_ptr<PSM> psm;
  } * ctx(sem);
  auto& t = *ctx;
  auto& fcu = ctx->fcu;
  auto& fca = ctx->fca;
  auto& fcn = ctx->fcn;
  auto& fci = ctx->fci;
  auto& fck = ctx->fck;
  auto& par = ctx->par;
  auto& psm = ctx->psm;

  if (sem("init")) {
    fcu.Reinit(m, 0);
    fca.Reinit(m, 0);
    fcn.Reinit(m, Vect(0.5));
    fci.Reinit(m, true);
  }
  if (sem.Nested()) {
    Plic plic{GRange<size_t>(1), &fcu,    &fca,    &fcn,  &fci,
              nullptr,           nullptr, nullptr, t.mebc};
    psm = UCurv<M>::CalcCurvPart(plic, par, &fck, m, m);
  }
  if (sem.Nested()) {
    psm->DumpParticles(&fca, &fcn, 0, 0);
  }
  if (sem("dump")) {
    m.Dump(&fck, "k");
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

set int px 1
set int py 1
set int pz 1
)EOF";

  MpiWrapper mpi(&argc, &argv);
  return RunMpiBasicString<M>(mpi, Run, conf);
}
