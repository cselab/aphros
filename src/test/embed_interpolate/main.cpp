// Created by Petr Karnakov on 05.07.2020
// Copyright 2020 ETH Zurich

#undef NDEBUG
#include <mpi.h>
#include <cassert>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

#include "distr/distrsolver.h"
#include "kernel/kernelmeshpar.h"
#include "parse/vars.h"

#include "dump/hdf.h"
#include "solver/embed.h"

template <class M_>
struct GPar {};

template <class M_>
class EmbedInterpolate : public KernelMeshPar<M_, GPar<M_>> {
 public:
  using P = KernelMeshPar<M_, GPar<M_>>;
  using M = M_;
  using Mesh = M;
  using Par = typename P::Par;
  static constexpr size_t dim = M::dim;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Sem = typename Mesh::Sem;

  EmbedInterpolate(Vars& var, const MyBlockInfo& block, Par& par);
  void Run() override;

 protected:
  using P::IsLead;
  using P::IsRoot;
  using P::m;
  using P::par_;
  using P::var;
};

template <class M>
EmbedInterpolate<M>::EmbedInterpolate(Vars& var, const MyBlockInfo& block, Par& par)
    : KernelMeshPar<M, Par>(var, block, par) {}

template <class M>
void EmbedInterpolate<M>::Run() {
  auto sem = m.GetSem("advection");
  if (sem()) {}
}

void Main(MPI_Comm comm, Vars& var) {
  using M = MeshStructured<double, 3>;
  using K = EmbedInterpolate<M>;
  using Par = typename K::Par;
  Par par;

  DistrSolver<M, K> ds(comm, var, par);
  ds.Run();
}

int main(int argc, const char** argv) {
  return RunMpi(argc, argv, Main);
}
