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


int main() {
  test_vars::Test();
}

/*
#include "Hydro.h"
void Main(MPI_Comm comm, bool loc) {
  // read config files, parse arguments, maybe init global fields
  using M = geom::geom3d::MeshStructured<Scal>;
  using K = Hydro<M>;
  using KF = HydroFactory<M>;
  using D = Distr;
  //using dx = D::Idx;

  KF kf;

  // Kernels have to know about Distr if they want to create Stage objects
  // However, Stage can be independent on Distr.
  // Comm() should put the list of fields to exchange somewhere
  // so that Distr could do communication.
  // Comm() must be independent on implementation of Distr.
  
  Idx b{1, 2, 2}; // number of blocks 
  Idx p{2, 1, 1}; // number of ranks
  const int es = 8;
  const int h = 1;
  const int bs = 16;
  
  // Initialize buffer mesh and make Hydro for each block.
  std::unique_ptr<Distr> d;
  if (loc) {
    d = CreateLocal(comm, kf, bs, b, p, es, h);
  } else {
    d = CreateCubism(comm, kf, bs, b, p, es, h);
  }

  while (!d->IsDone()) {
    // At each step, first exchange halos,
    // then call Kernel::operator() for each block
    d->Step();
  }
}

int main (int argc, const char ** argv) {
  int prov;
  MPI_Init_thread(&argc, (char ***)&argv, MPI_THREAD_MULTIPLE, &prov);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Comm comm;
  bool loc = true;
  if (loc) {
    MPI_Comm_split(MPI_COMM_WORLD, rank, rank, &comm);
    if (rank == 0) {
      Main(comm, loc);
    }
  } else {
    comm = MPI_COMM_WORLD;
    Main(comm, loc);
  }


  MPI_Finalize();	
  return 0;
}
*/
