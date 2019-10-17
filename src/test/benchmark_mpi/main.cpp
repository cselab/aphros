#undef NDEBUG
#include <iostream>

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
  for (auto i = 0; i < 10; ++i) {
    if (sem("std")) {
      s.v = {Scal(m.GetId())};
      using T = typename M::template OpCatT<Scal>;
      m.Reduce(std::make_shared<T>(&s.v));
    }
  }

  for (auto i = 0; i < 10; ++i) {
    if (sem("mpi")) {
      auto comm = MPI_COMM_WORLD;
      if (m.IsLead()) {
        auto& r = s.r;
        r = {0};
        int s = r.size(); // size local
        if (m.IsRoot()) {
          int sc; // size of communicator
          MPI_Comm_size(comm, &sc);

          std::vector<int> ss(sc); // size of r on all ranks

          MPI_Gather(&s, 1, MPI_INT, 
                     ss.data(), 1, MPI_INT, 
                     0, comm);

          int sa = 0; // size all
          std::vector<int> oo = {0}; // offsets
          for (auto& q : ss) {
            sa += q;
            oo.push_back(oo.back() + q);
          }
          oo.pop_back();
          assert(ss.size() == oo.size());

          std::vector<char> ra(sa); // result all

          MPI_Gatherv(r.data(), r.size(), MPI_CHAR,
                      ra.data(), ss.data(), oo.data(), MPI_CHAR,
                      0, comm);
        } else {
          MPI_Gather(&s, 1, MPI_INT, nullptr, 0, MPI_INT, 0, comm);
          MPI_Gatherv(r.data(), r.size(), MPI_CHAR,
                      nullptr, nullptr, nullptr, MPI_CHAR,
                      0, comm);
        }
      }
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

