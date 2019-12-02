#pragma once

#include <memory>
#include <mpi.h>
#include <stdexcept>
#include <functional>

#include "parse/vars.h"
#include "kernel/kernelmeshpar.h"
#include "distr.h"
#include "cubismnc.h"
#include "local.h"
#include "parse/parser.h"
#include "geom/vect.h"
#include "geom/block.h"
#include "distr/distrsolver.h"


std::string GetDefaultConf();

int RunMpi0(int argc, const char ** argv,
            std::function<void(MPI_Comm, Vars&)> r, std::istream& conf);

template <class M, class State, class Func>
int RunMpiBasic(int argc, const char** argv, Func func, std::string addconf) {
  struct Par { Func func; };
  class Basic : public KernelMeshPar<M, Par> {
   public:
    using P = KernelMeshPar<M, Par>;
    using P::par_;
    using P::m;
    using P::P;

    State s;
    void Run() {
      par_.func(m, s, this->var_mutable);
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
