// Created by Petr Karnakov on 01.10.2020
// Copyright 2020 ETH Zurich

#include <iostream>

#include "distr/distrbasic.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using MIdx = typename M::MIdx;

void Run(M& m, Vars&) {
  auto sem = m.GetSem();
  struct {
    FieldCell<Scal> fc_rhs;
    FieldCell<Scal> fcu;
  } * ctx(sem);
  auto& t = *ctx;
  if (sem()) {
    t.fc_rhs.Reinit(m, 0);
    for (auto c : m.CellsM()) {
      t.fc_rhs[c] = c.center[0];
    }
    t.fcu.Reinit(m, 0);
  }
  if (sem()) {
  }
}

int main(int argc, const char** argv) {
  std::string conf = R"EOF(
set int bx 2
set int by 2
set int bz 1

set int bsx 8
set int bsy 8
set int bsz 8

set int px 2
set int py 1
set int pz 1
)EOF";

  return RunMpiBasic<M>(argc, argv, Run, conf);
}
