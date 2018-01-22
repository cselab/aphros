/*
 *  main.cpp
 *  MPCFcluster
 *
 *  Created by Diego Rossinelli on 11/25/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <map>
#include <iomanip>
#include <iostream>
#include <mpi.h>

#include "Test_IO_Compressed.h"

using namespace std;

Simulation * sim = NULL;

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	const bool isroot = (myrank == 0);
	
	ArgumentParser parser(argc, (const char **)argv);	

	parser.set_strict_mode();
	
	if (!isroot)
		parser.mute();
    
	if( parser("-sim").asString() == "io" )
		sim = new Test_IO_Compressed(isroot, argc, (const char **)argv);
	else {
		if (isroot) printf("-sim value not recognized. Aborting.\n");
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	sim->setup();
	sim->dispose();
	
	delete sim;
	
	if (isroot) printf("Finishing...");
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	
	return 0;
}
