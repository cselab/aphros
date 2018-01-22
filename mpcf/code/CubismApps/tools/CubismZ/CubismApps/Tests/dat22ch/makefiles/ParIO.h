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
	MPI_Request request;	// for asynchronous writes

	string filename;
	int chunksize;
	
	MPI_ParIO()
	{
	}

	void Init(string name, size_t chunkbytes)
	{
		filename = name;
		chunksize = chunkbytes;
		
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
