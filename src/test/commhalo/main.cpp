#undef NDEBUG
#include <iostream>
#include <fstream>

#include "distr/distrbasic.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;

struct State {
  FieldCell<Scal> fc;
  std::ofstream out;
};

void Run(M& m, State& s, Vars&) {
  auto sem = m.GetSem();
  auto& fc = s.fc;
  auto& out = s.out;
  if (sem()) {
    out.open("o_" + std::to_string(m.GetId()) + ".log");
    fc.Reinit(m, -1);  // allocate memory for field fc
    out << "before id=" << m.GetId() << std::endl;
    for (auto c : m.AllCells()) {
      fc[c] = m.GetId();
      out << fc[c] << " ";
    }
    out << std::endl;
    m.Comm(&fc);      // exchange halo cells
  }
  if (sem()) {
    out << "after id=" << m.GetId() << std::endl;
    for (auto c : m.AllCells()) {
      fc[c] = m.GetId();
      out << fc[c] << " ";
    }
    out << std::endl;
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
set int bsz 1
)EOF";

  return RunMpiBasic<M, State>(argc, argv, Run, conf);
}

