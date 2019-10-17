#undef NDEBUG
#include <iostream>

#include "distr/distrbasic.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;

struct State {
  std::vector<Scal> v;
};

void Run(M& m, State& s, Vars& var) {
  (void) var;
  auto sem = m.GetSem();
  auto& v = s.v;
  for (auto i = 0; i < 10; ++i) {
    if (sem()) {
      v = {Scal(m.GetId())};
      using T = typename M::template OpCatT<Scal>;
      m.Reduce(std::make_shared<T>(&v));
    }
  }
  if (sem()) {
    if (m.IsRoot()) {
      m.TimerReport("timer.log");
    }
  }
  if (sem()) {}
}

int main(int argc, const char** argv) {
  std::string conf = R"EOF(
set int bx 2
set int by 2
set int bz 1

set int bsx 16
set int bsy 16
set int bsz 16

set int px 6
set int py 4
set int pz 8
  )EOF";

  return RunMpiBasic<M, State>(argc, argv, Run, conf);
}

