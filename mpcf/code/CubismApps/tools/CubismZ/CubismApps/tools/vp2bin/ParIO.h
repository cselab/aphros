#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <string>
using namespace std;

#pragma once

struct MPI_ParIO
{
	MPI_File outfile;
	MPI_Request request;	// for asychronous writes

	string filename;
	int chunksize;
	
	MPI_ParIO()
	{
	}

	void Init(string name, size_t nbytes)
	{
		filename = name;
		chunksize = nbytes;
		
		/* open file */
		int mode = MPI_MODE_CREATE | MPI_MODE_WRONLY;
		MPI_File_open( MPI_COMM_SELF, (char *)filename.c_str(), mode, MPI_INFO_NULL, &outfile );
	}

	void Write(char *data, int chunkid)
	{
		MPI_Status status;
		
		/* set file view */
		size_t offset = chunkid * chunksize; 
			
		MPI_File_write_at(outfile, offset, data, chunksize, MPI_BYTE, &status);
	}

	void Finalize()
	{
		MPI_File_close(&outfile);
	}

};


#if 0
#include <mpi.h>
#include "ParIO.h" 

int main(int argc, char *argv[])
{
	int rank, size;

	MPI_Init(&argc,&argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size); 

	/* Extract the original group handle */ 
	/* Divide tasks into two distinct groups based upon rank */ 

	int groupsize = 4;
	int reportfreq = 2;
	int layout = 0;

	MPI_ParIO pio;

	pio.Init((char *)"output.bin", groupsize, reportfreq, layout);	// 4, 2, 0

	int step;
	
	for (step = 0; step < 10; step++) {
		if (step% 2 == 0 && step > 0) {
			pio.Consolidate(step);
		}

//		float number = rank*10.0 + rand()%10  + step*100;
//		float number = rank;
		float number = step;
		pio.Notify(step, number);
	}

	pio.Finalize();

	MPI_Finalize();
	
	return 0;
} 
#endif
