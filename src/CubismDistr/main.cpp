/*
 *  main.cpp
 *  MPCFcluster
 *
 *  Created by Diego Rossinelli on 11/25/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <iostream>
#include <string>
#include <mpi.h>
#include <cassert>

#include "Vars.h"

#include "Hydro.h"
#include "Interp.h"

#include "../hydro/suspender.h"

void Main(MPI_Comm comm, bool loc, Vars& par) {
  // read config files, parse arguments, maybe init global fields
  using M = geom::geom3d::MeshStructured<Scal>;
  using K = Hydro<M>;
  using KF = HydroFactory<M>;
  using D = Distr;
  //using Idx = D::Idx;
  
  KF kf;

  const int es = 8;
  const int hl = par.Int["hl"];
  const int bs = 16;
  
  // Initialize buffer mesh and make Hydro for each block.
  std::unique_ptr<Distr> d;
  if (loc) {
    d = CreateLocal(comm, kf, bs, es, hl, par);
  } else {
    d = CreateCubism(comm, kf, bs, es, hl, par);
  }

  while (!d->IsDone()) {
    d->Step();
  }
}


int main (int argc, const char ** argv) {
  int prov;
  MPI_Init_thread(&argc, (char ***)&argv, MPI_THREAD_MULTIPLE, &prov);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Comm comm;

  Vars par;
  Interp ip(par);

  std::ifstream f("a.conf");
  ip.RunAll(f);
  ip.PrintAll();

  bool loc = par.Int["loc"];

  if (loc) {
    MPI_Comm_split(MPI_COMM_WORLD, rank, rank, &comm);
    if (rank == 0) {
      Main(comm, loc, par);
    }
  } else {
    comm = MPI_COMM_WORLD;
    Main(comm, loc, par);
  }


  MPI_Finalize();	
  return 0;
}
