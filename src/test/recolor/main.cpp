#undef NDEBUG
#include <iostream>

#include "distr/distrbasic.h"

using M = MeshStructured<double, 3>;

struct State {
};

void Run(M& m, State& s, Vars& var) {
  auto sem = m.GetSem();
  if (sem("init")) {
    (void) s;
    (void) var;
    std::cout << m.GetGlobalLength() << std::endl;
  }
}

int main(int argc, const char** argv) {
  return RunMpiBasic<M, State>(argc, argv, Run, "");
}

