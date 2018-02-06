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

/*

#include <ArgumentParser.h>
#include <Timer.h>

#include "Tests.h"

#define SIMP 1
#define HYDRO 2
#define CASE HYDRO

#if CASE == HYDRO
  #include "Test_Hydro.h"
#elif CASE == SIMP
  #include "Test_Simple.h"
  #include "SimpleStep.h"
#else
  #include "Test_SteadyStateMPI.h"
  #include "FlowStep_LSRK3MPI.h"
  typedef FlowStep_LSRK3MPI<GridMPI_t> TFlowStep;
#endif

using namespace std;

Simulation * sim = NULL;

int main (int argc, const char ** argv)
{
  int provided;
  MPI_Init_thread(&argc, (char ***)&argv, MPI_THREAD_MULTIPLE, &provided);

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  const bool isroot = (0 == myrank);

  ArgumentParser parser(argc, argv);

  if (parser.exist("simConf")) {
    parser.readFile(parser("simConf").asString());
  }

  parser.unset_strict_mode();

  if (!isroot)
    parser.mute();

  parser.set_strict_mode();

#if CASE == HYDRO
  sim = new Test_Hydro(MPI_COMM_WORLD);
#elif CASE == SIMP
  sim = new Test_Simple<GridMPI_t,SimpleStep<GridMPI_t>>(MPI_COMM_WORLD);
#else
  sim = new Test_SteadyStateMPI<GridMPI_t,FlowStep_LSRK3MPI<GridMPI_t>>(
      MPI_COMM_WORLD, parser);
#endif

  sim->setup();
  sim->run();
  delete sim;
  sim = NULL;

  MPI_Finalize();	

  return 0;
}

*/

#include <memory>
#include "Test_Hydro.h"

int main() {
  TestHydro sim(MPI_COMM_WORLD);
  sim.setup();
  sim.run();
}
