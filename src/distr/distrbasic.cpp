// Created by Petr Karnakov on 14.10.2019
// Copyright 2019 ETH Zurich

#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>

#include "distrsolver.h"
#include "util/distr.h"
#include "util/format.h"
#include "util/git.h"
#include "util/macros.h"
#include "util/subcomm.h"

std::string GetDefaultConf() {
#if USEFLAG(BACKEND_CUBISM)
  FORCE_LINK(distr_cubismnc);
#endif
#if USEFLAG(BACKEND_LOCAL)
  FORCE_LINK(distr_local);
#endif
#if USEFLAG(BACKEND_NATIVE)
  FORCE_LINK(distr_native);
#endif

#if USEFLAG(BACKEND_CUBISM)
  const std::string backend = "cubismnc";
#elif USEFLAG(BACKEND_NATIVE)
  const std::string backend = "native";
#else
  const std::string backend = "local";
#endif
  return util::Format(
      R"EOF(
set int px 1
set int py 1
set int pz 1
set int pw 1

set int bx 1
set int by 1
set int bz 1
set int bw 1

set int bsx 16
set int bsy 16
set int bsz 16
set int bsw 16

set int CHECKNAN 1
set int dim 3

set string dumpformat default

set int comm_size 8

set string backend {}
set int histogram 0
set int openmp 0
set int mpi_compress_msg 0

set int max_step 1
set int num_frames 1

set int hl 2
set int hypre_print 0
set double hypre_symm_tol 1e-12
set double hypre_vort_tol 1e-12
set double hypre_gen_tol 1e-12
set int periodic 1
set int hypre_periodic_x 1
set int hypre_periodic_y 1
set int hypre_periodic_z 1
set string hypre_gen_solver gmres

set int linsolver_gen_maxnorm 0
set int linsolver_symm_maxnorm 0
set int linsolver_vort_maxnorm 0

set int loc_periodic_x 1
set int loc_periodic_y 1
set int loc_periodic_z 1
set int loc_periodic_w 1
set int loc_maxcomm 16

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
)EOF",
      backend);
}

int RunMpiKernel(
    MpiWrapper& mpi, std::function<void(MPI_Comm, Vars&)> kernel,
    std::istream& conf) {
  Vars var;
  Parser(var).ParseStream(conf);

  const std::string backend = var.String["backend"];
  if (backend == "local") {
    fassert_equal(
        mpi.GetCommSize(), 1, "\nBackend 'local' requires a single rank.\n");
  }

  try {
    kernel(mpi.GetComm(), var);
  } catch (const std::exception& e) {
    std::cerr << FILELINE + "\nabort after throwing exception\n"
              << e.what() << '\n';
    std::terminate();
  }

  return 0;
}
