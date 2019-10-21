#undef NDEBUG
#include <iostream>

#include "distr/distrbasic.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;

struct State {
  FieldCell<Scal> fc;
};

void Run(M& m, State& s, Vars&) {
  auto sem = m.GetSem();
  if (sem()) {
    fc.Reinit(m, 0);  // allocate memory for field fc
    std::cout << "before id=" << m.GetId() << std::endl;
    for (auto c : m.Cells()) {
      fc[c] = m.GetId();
      std::cout << fc[c] << " ";
    }
    std::cout << std::endl;
    m.Comm(&fc);      // exchange halo cells
  }
  if (sem()) {
    std::cout << "after id=" << m.GetId() << std::endl;
    for (auto c : m.Cells()) {
      fc[c] = m.GetId();
      std::cout << fc[c] << " ";
    }
    std::cout << std::endl;
  }
}

int main(int argc, const char** argv) {
  std::string conf = R"EOF(
set int px 2
set int py 1
set int pz 1

set int bx 2
set int by 1
set int bz 1

set int bsx 8
set int bsy 8
set int bsz 8
)EOF";

  return RunMpiBasic<M, State>(argc, argv, Run, conf);
}

