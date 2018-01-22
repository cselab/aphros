/*
 *  main.cpp
 *  MPCFcluster
 *
 *  Created by Panagiotis Hadjidoukas on 4/4/16.
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */

#include <iostream>
#include <mpi.h>

//#include "Test_IO_Compressed.h"
#include "Test_Statistics.h"
#include "ArgumentParser.h"

using namespace std;

Simulation * sim = NULL;

int main (int argc, char ** argv)
{
	MPI_Init(&argc, &argv);
	
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	const bool isroot = myrank == 0;
	
	ArgumentParser parser (argc, (const char **)argv);
  
	parser.set_strict_mode();
	
	if (!isroot)
		parser.mute();
  
	if ( parser("-sim").asString() == "statistics" )
		sim = new Test_Statistics (isroot, argc, (const char **)argv);
	//if ( parser("-sim").asString() == "io" )
	//	sim = new Test_IO_Compressed(isroot, argc, (const char **)argv);
	else {
		if (isroot)
		{
			printf("-sim value not recognized. Aborting.\n");
			printf("usage: %s -inputpath <path> -outputpath <path> -stepid <number> [-swap] [-restart] [-wtype_read <type1>] [-wtype_write <type2>]\n", argv[0]);
		}
		MPI_Abort(MPI_COMM_WORLD, 1);
	} 
 
	sim->setup();
	
	sim->dispose();
	
	delete sim;
	
	sim = NULL;
	
	if (isroot) printf("Finishing...");
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	
	return 0;
}
