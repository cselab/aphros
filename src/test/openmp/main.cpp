// Created by Petr Karnakov on 25.11.2019
// Copyright 2019 ETH Zurich

#undef NDEBUG

#ifdef _OPENMP
#include <omp.h>
#endif

#include <sched.h>
#include <iostream>

#include "distr/distrbasic.h"

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;
using MIdx = typename M::MIdx;

void Run(M& m, Vars&) {
#ifdef _OPENMP
  const static double wtime0 = omp_get_wtime();
  auto sem = m.GetSem();
  auto blockid = [&m]() -> MIdx {
    auto& bc = m.GetInBlockCells();
    return bc.GetBegin() / bc.GetSize();
  };
  if (sem()) {
    if (m.IsRoot()) {
      std::cerr << "\nRun()" << std::endl;
    }
  }
  if (sem()) {
    volatile double a = 10.;
    for (size_t i = 0; i < (1 << 28); ++i) {
      a = std::sqrt(a);
    }
#pragma omp critical
    {
      std::cerr << "block=" << std::setw(3) << blockid();
      std::cerr << std::setw(10) << " thread=" << std::setw(2)
                << omp_get_thread_num();
      std::cerr << std::setw(8) << " cpu=" << std::setw(2) << sched_getcpu();
      std::cerr << std::setw(8) << " t=" << std::setw(8)
                << omp_get_wtime() - wtime0;
      std::cerr << std::endl;
    }
  }
#endif
}

int main(int argc, const char** argv) {
  std::string conf = R"EOF(
set int bx 4
set int by 1
set int bz 1

set int bsx 16
set int bsy 16
set int bsz 16

set int px 2
set int py 1
set int pz 1

set string backend cubismnc
set int loc_maxcomm 16
set int verbose_openmp 1
)EOF";

  MpiWrapper mpi(&argc, &argv);
  return RunMpiBasicString<M>(mpi, Run, conf);
}
