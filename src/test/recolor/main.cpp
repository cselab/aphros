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

int RunMpi0(int argc, const char ** argv,
           std::function<void(MPI_Comm, Vars&)> r) {
  int prov;
  MPI_Init_thread(&argc, (char ***)&argv, MPI_THREAD_MULTIPLE, &prov);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  Vars var;   // parameter storage
  Parser ip(var); // parser

  std::stringstream f;  // config

  f << R"foo(
set int px 1
set int py 1
set int pz 1

set int bx 1
set int by 1
set int bz 1

set int bsx 16
set int bsy 16
set int bsz 16

set int CHECKNAN 1
set int dim 3

set int comm_size 8

set int loc 0

set int max_step 1
set int num_frames 1

set int hl 2
set int hypre_print 2
set double hypre_symm_tol 1e-12
set double hypre_vort_tol 1e-12
set double hypre_gen_tol 1e-12
set int periodic 1
set int hypre_periodic_x 1
set int hypre_periodic_y 1
set int hypre_periodic_z 1
set string hypre_gen_solver gmres

set int loc_periodic_x 1
set int loc_periodic_y 1
set int loc_periodic_z 1

set int verbose 0
set int output 0
set int verbose_stages 0
set int verbose_time 0

set int iter 1

set double extent 1.

set string hypre_symm_solver pcg
set string hypre_vort_solver pcg
set int hypre_symm_maxiter 100
set int hypre_vort_maxiter 100
set int hypre_gen_maxiter 30
)foo";

  ip.RunAll(f);

  bool loc = var.Int["loc"];

  MPI_Comm comm;
  if (loc) {
    MPI_Comm_split(MPI_COMM_WORLD, rank, rank, &comm);
    if (rank == 0) {
      r(comm, var);
    }
  } else {
    comm = MPI_COMM_WORLD;
    r(comm, var);
  }

  MPI_Finalize();
  return 0;
}

template <class M_, class State, class Func>
int RunMpiBasic(int argc, const char** argv, Func func) {
  struct Par { Func func; };
  class Basic : public KernelMeshPar<M_, Par> {
   public:
    using P = KernelMeshPar<M_, Par>;
    using P::var;
    using P::par_;
    using P::m;
    using P::P;

    State s;
    void Run() {
      par_.func(m, s, var);
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

  Main main(func);
  return RunMpi0(argc, argv, main);
}


int main(int argc, const char** argv) {
  return RunMpiBasic<M, State>(argc, argv, Run);
}

