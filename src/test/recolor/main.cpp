#undef NDEBUG
#include <iostream>
#include <string>
#include <mpi.h>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>
#include <cmath>

#include "geom/mesh.h"
#include "kernel/kernelmeshpar.h"
#include "distr/distrsolver.h"
#include "distr/distrbasic.h"
#include "util/suspender.h"
#include "solver/solver.h"
#include "linear/linear.h"
#include "solver/pois.h"

using M = MeshStructured<double, 3>;

struct State {
};

void Run(M& m, State& s, Vars& var) {
  auto sem = m.GetSem("Comm");
  if (sem("init")) {
    (void) s;
    (void) var;
    std::cout << m.GetGlobalLength() << std::endl;
  }
}

int main(int argc, const char** argv) {
  return RunMpiBasic<M, State>(argc, argv, Run, "");
}

