// Created by Petr Karnakov on 17.10.2019
// Copyright 2019 ETH Zurich

#undef NDEBUG
#include <iostream>

#include "distr/distrbasic.h"

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;

void Run(M& m, Vars&) {
  auto sem = m.GetSem();
  struct {
    std::vector<Scal> v;
    std::vector<char> r;
  } * ctx(sem);
  auto& v = ctx->v;
  auto& r = ctx->r;
  for (auto i = 0; i < 10; ++i) {
    if (sem("std")) {
      v = {Scal(m.GetId())};
      m.Reduce(&v, Reduction::concat);
    }
  }

  for (auto i = 0; i < 10; ++i) {
    if (sem("mpi")) {
      auto comm = MPI_COMM_WORLD;
      if (m.IsLead()) {
        r = {0};
        int s = r.size(); // size local
        if (m.IsRoot()) {
          int sc; // size of communicator
          MPI_Comm_size(comm, &sc);

          std::vector<int> ss(sc); // size of r on all ranks

          MPI_Gather(&s, 1, MPI_INT, ss.data(), 1, MPI_INT, 0, comm);

          int sa = 0; // size all
          std::vector<int> oo = {0}; // offsets
          for (auto& q : ss) {
            sa += q;
            oo.push_back(oo.back() + q);
          }
          oo.pop_back();
          assert(ss.size() == oo.size());

          std::vector<char> ra(sa); // result all

          MPI_Gatherv(
              r.data(), r.size(), MPI_CHAR, ra.data(), ss.data(), oo.data(),
              MPI_CHAR, 0, comm);
        } else {
          MPI_Gather(&s, 1, MPI_INT, nullptr, 0, MPI_INT, 0, comm);
          MPI_Gatherv(
              r.data(), r.size(), MPI_CHAR, nullptr, nullptr, nullptr, MPI_CHAR,
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
  if (sem()) {
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

set int px 6
set int py 4
set int pz 8
)EOF";

  MpiWrapper mpi(&argc, &argv);
  return RunMpiBasicString<M>(mpi, Run, conf);
}
