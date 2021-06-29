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
using MIdx = typename M::MIdx;

void Run(M& m, Vars&) {
  auto sem = m.GetSem();

  const GRange<size_t> layers{1};
  struct {
    FieldCell<Scal> fcu; // volume fraction
    FieldCell<Scal> fca; // plane constant
    FieldCell<Vect> fcn; // normal
    FieldCell<bool> fci; // interface mask
    FieldCell<Scal> fck; // curvature
    MapEmbed<BCondAdvection<Scal>> mebc;
    typename PartStrMeshM<M>::Par par;
    std::unique_ptr<curvature::Particles<M>> curv_part;
  } * ctx(sem);
  auto& t = *ctx;

  if (sem("init")) {
    t.fcu.Reinit(m, 0);
    t.fca.Reinit(m, 0);
    t.fcn.Reinit(m, Vect(0.5));
    t.fci.Reinit(m, true);
    t.curv_part.reset(new curvature::Particles<M>(m, t.par, layers));
  }
  if (sem.Nested()) {
    generic::Plic<Vect> plic{layers,  &t.fcu,  &t.fca,  &t.fcn, &t.fci,
                             nullptr, nullptr, nullptr, t.mebc};
    t.curv_part->CalcCurvature(&t.fck, plic, m, m);
  }
  if (sem.Nested()) {
    t.curv_part->GetParticles()->DumpParticles(&t.fca, &t.fcn, 0, 0);
  }
  if (sem("dump")) {
    m.Dump(&t.fck, "k");
  }
  if (sem()) {
  }
}

int main(int argc, const char** argv) {
  MpiWrapper mpi(&argc, &argv);
  const MIdx mesh_size(16);
  Subdomains<MIdx> sub(mesh_size, mesh_size, mpi.GetCommSize());
  std::string conf = sub.GetConfig();
  return RunMpiBasicString<M>(mpi, Run, conf);
}
