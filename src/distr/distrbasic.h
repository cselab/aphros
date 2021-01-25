// Created by Petr Karnakov on 14.10.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <mpi.h>
#include <fstream>
#include <functional>
#include <memory>
#include <sstream>
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
#include "util/distr.h"

std::string GetDefaultConf();

int RunMpiKernel(
    MpiWrapper& mpi, std::function<void(MPI_Comm, Vars&)> r,
    std::istream& conf);

template <class M>
int RunMpiBasicString(
    MpiWrapper& mpi, std::function<void(M& m, Vars&)> kernel,
    std::string addconf) {
  using Kernel = std::function<void(M & m, Vars&)>;
  struct Par {
    Kernel kernel;
  };
  class Basic : public KernelMeshPar<M, Par> {
   public:
    using P = KernelMeshPar<M, Par>;
    using P::m;
    using P::P;
    using P::par_;

    void Run() {
      par_.kernel(m, this->var_mutable);
    }
  };

  struct Main {
    Main(Kernel kernel_) : kernel(kernel_) {}
    void operator()(MPI_Comm comm, Vars& var) {
      Par par{kernel};
      DistrSolver<M, Basic> ds(comm, var, par);
      ds.Run();
    }
    Kernel kernel;
  };

  std::stringstream conf;
  conf << GetDefaultConf();
  conf << "\n" << addconf;
  Main main(kernel);
  return RunMpiKernel(mpi, main, conf);
}

template <class M>
int RunMpiBasicFile(
    MpiWrapper& mpi, std::function<void(M& m, Vars&)> kernel,
    std::string confpath = "a.conf") {
  std::ifstream file(confpath);
  if (!file.good()) {
    throw std::runtime_error(
        FILELINE + ": Can't open config file '" + confpath + "'");
  }
  std::stringstream buf;
  buf << file.rdbuf();
  return RunMpiBasicString<M>(mpi, kernel, buf.str());
}
