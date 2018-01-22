/* File:   main.cpp */
/* Date:   July 2015 */
/* Author: Ursula Rasthofer */
/* Tag:    spectral space evaluations of cubism results */
/* Copyright 2015 ETH Zurich. All Rights Reserved. */

#include <iostream>
#include <mpi.h>

#include "ArgumentParser.h"
#include "evalBase.h"

#include "eval_hit/evalHit.h"

using namespace std;

int main(int argc, char *argv[])
{
  // initialize MPI
  MPI_Init(&argc, &argv);
  // get rank
  const int myrank = MPI::COMM_WORLD.Get_rank();

  // read configuration
  ArgumentParser myparser(argc, (const char **)argv);
  myparser.readFile(myparser("-configFile").asString("./config.dat"));
  // and print arguments
  if (myrank == 0)
    myparser.print_args();

   // evaluate in spectral space
  evalBase * evaluation = 0;

  // select case
  if (myparser("sim").asString() == "hit" )
    evaluation = new evalHit(myparser);
  else
  {
    if (myrank == 0)
      printf("sim value not recognized. Aborting.\n");
    abort();
  }

  // do evaluations in spectral space
  evaluation->compute();

  if (myrank == 0)
    cout << "FINISHED" << endl;

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}
