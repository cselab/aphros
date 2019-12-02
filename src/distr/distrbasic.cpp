#include "distrsolver.h"
#include "util/git.h"

std::string GetDefaultConf() {
  return R"foo(
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

set string dumpformat default

set int comm_size 8

set string backend cubismnc
set int histogram 0

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
set int verbose_openmp 0

set int iter 1

set double extent 1.

set string hypre_symm_solver pcg
set string hypre_vort_solver pcg
set int hypre_symm_maxiter 100
set int hypre_vort_maxiter 100
set int hypre_gen_maxiter 30
)foo";
}

int RunMpi0(int argc, const char ** argv,
            std::function<void(MPI_Comm, Vars&)> r, std::istream& conf) {
  int prov;
  MPI_Init_thread(&argc, (char ***)&argv, MPI_THREAD_MULTIPLE, &prov);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  Vars var;   // parameter storage
  Parser ip(var); // parser

  ip.RunAll(conf);

  std::string be = var.String["backend"];

  MPI_Comm comm;
  if (be == "local") {
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



