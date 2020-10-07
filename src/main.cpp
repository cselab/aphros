// Created by Petr Karnakov on 30.05.2018
// Copyright 2018 ETH Zurich

#include "distr/distrsolver.h"
#include "kernel/hydro.h"
#include "kernel/kernelmeshpar.h"

void Main(MPI_Comm comm, Vars& var) {
  FORCE_LINK(init_contang);
  FORCE_LINK(init_vel);

  using M = MeshStructured<double, 3>;
  using K = Hydro<M>;
  using Par = typename K::Par;
  Par par;

  DistrSolver<M, K> ds(comm, var, par);
  ds.Run();
}

int main(int argc, const char** argv) {
  return RunMpi(argc, argv, Main);
}
