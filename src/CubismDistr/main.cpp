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

#include "../hydro/suspender.h"
#include "Vars.h"
#include "Hydro.h"
#include "Interp.h"
#include "ICubism.h"
#include "ILocal.h"

void Main(MPI_Comm comm, bool loc, Vars& par) {
  using Scal = double;
  using M = geom::MeshStructured<Scal, 3>;
  using KF = HydroFactory<M>;
  
  KF kf;

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  bool isroot = (!rank);

  std::unique_ptr<Distr> d;
  if (loc) {
    std::cerr << "Create Local on " << size << " ranks" << std::endl;
    d.reset(CreateLocal(comm, kf, par));
  } else {
    std::cerr << "Create Cubism on " << size << " ranks" << std::endl;
    d.reset(CreateCubism(comm, kf, par));
  }

  std::cerr << "Run main loop" << std::endl;
  d->Run();
}


int main (int argc, const char ** argv) {
  int prov;
  MPI_Init_thread(&argc, (char ***)&argv, MPI_THREAD_MULTIPLE, &prov);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  bool isroot = (!rank);

  Vars par;   // parameter storage
  Interp ip(par); // interpretor (parser)

  std::string fn = "a.conf";
  if (argc == 1) {
    // nop
  } else if (argc == 2) {
    fn = argv[1];
  } else {
    if (isroot) {
      std::cerr << "usage: " << argv[0] << " [a.conf]" << std::endl;
    }
    return 1;
  }

  if (isroot) {
    std::cerr << "Loading config from '" << fn << "'" << std::endl;
  }

  std::ifstream f(fn);  // config file
  // Read file and run all commands
  ip.RunAll(f);   

  // Print vars on root
  if (isroot) {            
    std::cerr << "\n=== config begin ===" << std::endl;
    ip.PrintAll(std::cerr);
    std::cerr << "=== config end ===\n" << std::endl;
  }

  bool loc = par.Int["loc"];

  MPI_Comm comm;
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
