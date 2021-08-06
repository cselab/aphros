// Created by Sergey Litvinov on 01.03.2021
// Copyright 2021 ETH Zurich

#include "aphros_c.h"
#include "distr/distrsolver.h"
#include "kernel/hydro.h"
#include "kernel/kernelmeshpar.h"
#include "overlap/overlap.h"

template <size_t dim>
static void Run(MPI_Comm comm, Vars& var) {
  using M = MeshCartesian<double, dim>;
  typename Hydro<M>::Par par;

  DistrSolver<M, Hydro<M>> ds(comm, var, par);
  ds.Run();
}

static void Main(MPI_Comm comm, Vars& var) {
  FORCE_LINK(init_contang);
  FORCE_LINK(init_vel);

  const int dim = var.Int("spacedim", 3);
  switch (dim) {
#if USEFLAG(DIM1)
    case 1:
      Run<1>(comm, var);
      break;
#endif
#if USEFLAG(DIM2)
    case 2:
      Run<2>(comm, var);
      break;
#endif
#if USEFLAG(DIM3)
    case 3:
      Run<3>(comm, var);
      break;
#endif
#if USEFLAG(DIM4)
    case 4:
      Run<4>(comm, var);
      break;
#endif
    default:
      fassert(false, "Unknown dim=" + std::to_string(dim));
  }
}

int aphros_Main(int argc, const char** argv) {
  return RunMpi(argc, argv, Main);
}

int aphros_GetSphereOverlap(const double *x, const double *h, const double *c, double s,
			    /**/ double *ans) {
  *ans = GetSphereOverlap(generic::Vect<double, 3>(x),
			  generic::Vect<double, 3>(h),
			  generic::Vect<double, 3>(c),
			  s);
  return 0;
}
