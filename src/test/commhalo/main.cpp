// Created by Petr Karnakov on 21.10.2019
// Copyright 2019 ETH Zurich

#undef NDEBUG
#include <fstream>
#include <iostream>

#include "distr/distrbasic.h"

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;

void Run(M& m, Vars&) {
  auto sem = m.GetSem();
  struct {
    FieldCell<Scal> fc;
    std::ofstream out;
  } * ctx(sem);
  auto& fc = ctx->fc;
  auto& out = ctx->out;
  auto& ic = m.GetIndexCells();
  auto& bc = m.GetAllBlockCells();
  if (sem("local")) {
    out.open("o_" + std::to_string(m.GetId()) + ".log");
    fc.Reinit(m, -1); // allocate memory for field fc
    out << "before id=" << m.GetId() << std::endl;
    for (auto c : m.Cells()) { // traverse internal cells
      fc[c] = m.GetId();
    }
    for (auto w : bc) { // traverse all cells with multi-index
      IdxCell c = ic.GetIdx(w);
      out << w << " " << fc[c] << std::endl;
    }
    out << std::endl;
    m.Comm(&fc); // exchange halo cells
  }
  if (sem("comm")) {
    out << "after id=" << m.GetId() << std::endl;
    for (auto w : bc) {
      IdxCell c = ic.GetIdx(w);
      out << w << " " << fc[c] << std::endl;
    }
    out << std::endl;
  }
  if (sem()) {
  }
}

int main(int argc, const char** argv) {
  std::string conf = R"EOF(
# ranks
set int px 2
set int py 1
set int pz 1

# blocks per rank
set int bx 2
set int by 1
set int bz 1

# cells per block
set int bsx 8
set int bsy 8
set int bsz 1

set int hl 2

set int verbose_time 1
set int verbose_stages 1
)EOF";

  MpiWrapper mpi(&argc, &argv);
  return RunMpiBasicString<M>(mpi, Run, conf);
}
