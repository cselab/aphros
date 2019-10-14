#undef NDEBUG
#include <iostream>

#include "distr/distrbasic.h"
#include "func/init_u.h"
#include "solver/solver.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;

struct State {
  FieldCell<Scal> fcu;
};

void Run(M& m, State& s, Vars& var) {
  (void) var;
  auto sem = m.GetSem();
  auto& fcu = s.fcu;
  if (sem()) {
    fcu.Reinit(m);
    auto init = CreateInitU<M>(var, m.IsRoot());
    init(fcu, m);
  }
  for (size_t i = 0; i < 10; ++i) {
    if (sem.Nested()) {
      solver::Smoothen(fcu, MapFace<std::shared_ptr<solver::CondFace>>(), m, 1);
    }
    if (sem()) {
      m.Dump(&fcu, "u");
    }
  }
  if (sem()) {}
}

int main(int argc, const char** argv) {
  std::string conf = R"EOF(
set int bx 2
set int by 2
set int bz 1

set int bsx 8
set int bsy 8
set int bsz 8

set string init_vf list
set int list_ls 1
set string list_path b.dat
  )EOF";

  return RunMpiBasic<M, State>(argc, argv, Run, conf);
}

