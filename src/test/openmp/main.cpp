#undef NDEBUG
#include <iostream>
#include <omp.h>
#include <sched.h>

#include "distr/distrbasic.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;

struct State {
  std::vector<Scal> v;
  std::vector<char> r;
  Scal a;
};

void Run(M& m, State& s, Vars& var) {
  (void) var;
  auto sem = m.GetSem();
  if (sem()) {
    volatile double a = 10.;
    for (size_t i = 0; i < (1 << 28); ++i) {
      a = std::sqrt(a);
    }
    std::cerr << "thread=" << omp_get_thread_num()
        << " cpu=" << sched_getcpu()
        << std::endl;
  }
}

int main(int argc, const char** argv) {
  std::string conf = R"EOF(
set int bx 2
set int by 2
set int bz 1

set int bsx 16
set int bsy 16
set int bsz 16

set int px 1
set int py 1
set int pz 1

set string backend cubismnc
)EOF";

  return RunMpiBasic<M, State>(argc, argv, Run, conf);
}

