// Created by Petr Karnakov on 14.10.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <mpi.h>
#include <functional>
#include <memory>
#include <stdexcept>

#include "cubismnc.h"
#include "distr.h"
#include "distr/distrsolver.h"
#include "geom/block.h"
#include "geom/vect.h"
#include "kernel/kernelmeshpar.h"
#include "local.h"
#include "parse/parser.h"
#include "parse/vars.h"

std::string GetDefaultConf();

int RunMpi0(
    int argc, const char** argv, std::function<void(MPI_Comm, Vars&)> r,
    std::istream& conf);

template <class M, class Func>
int RunMpiBasic(int argc, const char** argv, Func func, std::string addconf) {
  struct Par {
    Func func;
  };
  class Basic : public KernelMeshPar<M, Par> {
   public:
    using P = KernelMeshPar<M, Par>;
    using P::m;
    using P::P;
    using P::par_;

    void Run() {
      par_.func(m, this->var_mutable);
    }
  };

  struct Main {
    Main(Func func) : func(func) {}
    void operator()(MPI_Comm comm, Vars& var) {
      Par par{func};
      DistrSolver<M, Basic> ds(comm, var, par);
      ds.Run();
    }
    Func func;
  };

  std::stringstream conf;
  conf << GetDefaultConf();
  conf << "\n" << addconf;
  Main main(func);
  return RunMpi0(argc, argv, main, conf);
}
