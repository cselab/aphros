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

#include <ArgumentParser.h>
#include <Timer.h>

#include "Tests.h"
#include "Test_SteadyStateMPI.h"
//#include "FlowStep_LSRK3MPI.h"
#include "SimpleStep.h"

//typedef FlowStep_LSRK3MPI<GridMPI_t> TFlowStep;

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

  //sim = new Test_SteadyStateMPI<GridMPI_t,FlowStep_LSRK3MPI<GridMPI_t>>(
  //    MPI_COMM_WORLD, parser);

  sim = new Test_SteadyStateMPI<GridMPI_t,SimpleStep<GridMPI_t>>(
      MPI_COMM_WORLD, parser);

  sim->setup();
  sim->run();
  delete sim;
  sim = NULL;

  MPI_Finalize();	

  return 0;
}
