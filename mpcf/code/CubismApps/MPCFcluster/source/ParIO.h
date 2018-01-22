#pragma once

#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

#define MAX_REPORTFREQ	512

struct MPI_ParIO_Group
{
    MPI_Comm world_comm;

	MPI_Group orig_group;	// from world_comm
	int rank, size;

	MPI_Group new_group;	// groups for gathering
	MPI_Comm new_comm;
	int new_rank, new_size;

	MPI_Group master_group;	// group of masters
	MPI_Comm master_comm;
	int master_rank, master_size;

	int groupsize;
	int reportfreq;
	int layout;	// 0: chunked, 1: ordered, 2: chunked+async, 3: ordered+async
	MPI_Datatype filetype, buftype;	// for ordered

	int ngroups; // = size / groupsize;

	MPI_ParIO_Group()
	{
	}

	void Init(int _groupsize, int _reportfreq, int _layout, MPI_Comm comm)
	{
        world_comm = comm;

		MPI_Comm_group(world_comm, &orig_group);

		MPI_Comm_rank(world_comm, &rank);
		MPI_Comm_size(world_comm, &size);

		groupsize = _groupsize;
		reportfreq = _reportfreq;
		layout = _layout;

		int range[1][3];
//		int ngroups = size / groupsize;
//		printf("ngroups = %d\n", ngroups); fflush(0);

		int mygroup = rank / groupsize;
		int first = mygroup * groupsize;
		int last = first + groupsize-1;
		if (last >= size) last = size-1;	// correct groupsize for the last group

		range[0][0] = first;
		range[0][1] = last;
		range[0][2] = 1;

//		printf("[%d] first = %d last = %d\n", rank, first, last); fflush(0);
		MPI_Group_range_incl(orig_group, 1, range, &new_group);
		MPI_Group_rank (new_group, &new_rank);
		MPI_Group_size (new_group, &new_size);

		MPI_Comm_create(world_comm, new_group, &new_comm);

		int range_master[1][3];

		range_master[0][0] = 0;
		range_master[0][1] = size-1;
		range_master[0][2] = groupsize;

		master_comm = MPI_COMM_NULL;
		master_rank = -1;
		master_size = 0;

		MPI_Group_range_incl(orig_group, 1, range_master, &master_group);

		MPI_Comm_create(world_comm, master_group, &master_comm);
		MPI_Group_rank (master_group, &master_rank);
		MPI_Group_size (master_group, &master_size);

		/* ordered layout and single write call */
		if ((layout == 1) || (layout == 3)) {
			/* create and commit filetype */
			MPI_Type_vector(reportfreq, groupsize, size, MPI_FLOAT, &filetype);
			MPI_Type_commit(&filetype);

			/* create and commit buftype */
			MPI_Type_contiguous(groupsize * reportfreq, MPI_FLOAT, &buftype);
			MPI_Type_commit(&buftype);
		}

	}

};

struct MPI_ParIO
{

	MPI_File outfile;
	MPI_Request request;	// for asychronous writes

	string filename;

//	float numbers[MAX_REPORTFREQ];	// vector!
	float *numbers;
	int numid;	// vector

	int mygstep;
	long jobid;

	int async_counter;
	int lastcall;

	float *all_numbers;

	MPI_ParIO_Group *group;

	int groupsize, reportfreq, layout;
	int	rank, new_rank, master_rank;
	int	size, new_size, master_size;
	MPI_Datatype filetype, buftype;	// for ordered


	MPI_ParIO(int bsize=512): mygstep(0), async_counter(0), lastcall(0), all_numbers(NULL)
	{
//		printf("allocating %d intermediate buffer size\n", bsize);
		numbers = (float *)malloc(bsize*sizeof(float));
	}

	void Init(string name, MPI_ParIO_Group *_group)
	{
		group = _group;

#if 1	/* copy from group */
		groupsize = group->groupsize;
		reportfreq = group->reportfreq;
		layout = group->layout;

		rank = group->rank;
		new_rank = group->new_rank;
		master_rank = group->master_rank;

		size = group->size;
		new_size = group->new_size;
		master_size = group->master_size;

		filetype = group->filetype;
		buftype = group->buftype;

/*
		printf("%d %d %d \n %d %d %d \n %d %d %d\n",
			groupsize, reportfreq, layout,
			rank, new_rank, master_rank,
			size, new_size, master_size);
*/
#endif

		filename = name;

#if 0
		if (rank == 0) {
			jobid = getpid();
		}
		MPI_Bcast(&jobid, 1, MPI_LONG, 0, group->world_comm);

		std::stringstream ss;
		ss << "." << jobid;
		filename.append(ss.str());
#endif

#if 1
		if (rank == 0) {
			MPI_File testfile;

			int mode = MPI_MODE_RDONLY;
			int rc = MPI_File_open( MPI_COMM_SELF, (char *)filename.c_str(), mode, MPI_INFO_NULL, &testfile);
			if (!rc) {
 				printf("Warning: Hist file %s already exists and will be deleted (error code %d).\n", (char *)filename.c_str(), rc);fflush(stdout);
				rc = MPI_File_delete((char *)filename.c_str(), MPI_INFO_NULL);
				if (rc) {
					printf("Unable to delete file %s, error code %d\n", (char *)filename.c_str(), rc);fflush(stdout);
				}
				MPI_File_close(&testfile);
			}
		}
		MPI_Barrier(group->world_comm);
#endif


		all_numbers = (float *)malloc(groupsize*reportfreq*sizeof(float));	// explicit allocation of data buffer

		numid = 0;

		if (new_rank == 0) {	// parallel io
//			printf("I am here [master_comm = 0x%lx]!\n", master_comm);

			/* open file */
			int mode = MPI_MODE_CREATE | MPI_MODE_WRONLY;
			//MPI_File_open( master_comm, (char *)filename.c_str(), mode, MPI_INFO_NULL, &outfile );
			MPI_File_open( MPI_COMM_SELF, (char *)filename.c_str(), mode, MPI_INFO_NULL, &outfile );

		}
	}

	void Notify(float number)
	{
		numbers[numid] = number;
		numid++;
	}

	void _consolidate_chunked_async(int gstep)
	{
		async_counter++;

		int step = (gstep - reportfreq)/reportfreq;	// normalized step

		numid = 0;	// reset the vector ;-)


		if (new_rank == 0) {
			MPI_Status status;
			if (async_counter > 1) MPI_Wait(&request, &status);
		}

		MPI_Gather(&numbers[0], reportfreq, MPI_FLOAT, all_numbers, reportfreq, MPI_FLOAT, 0, group->new_comm);	// 0: master of the group (new_rank == 0)

		if (new_rank == 0) {	// master: perform parallel io
			MPI_Status status;

			size_t base = step * size * reportfreq;
			size_t offset = master_rank * groupsize * reportfreq;

			MPI_File_iwrite_at(outfile, (base + offset)*sizeof(float), &all_numbers[0], reportfreq*groupsize, MPI_FLOAT, &request);
			if (lastcall) {
				MPI_Wait(&request, &status);
			}
		}

	}

	void _consolidate_chunked(int gstep)
	{
		int step = (gstep - reportfreq)/reportfreq;	// normalized step

		numid = 0;	// reset the vector ;-)

		MPI_Gather(&numbers[0], reportfreq, MPI_FLOAT, all_numbers, reportfreq, MPI_FLOAT, 0, group->new_comm);	// 0: master of the group (new_rank == 0)

		if (new_rank == 0) {	// master: perform parallel io
			MPI_Status status;

			size_t base = step * size * reportfreq;
			size_t offset = master_rank * groupsize * reportfreq;

			MPI_File_write_at(outfile, (base + offset)*sizeof(float), &all_numbers[0], reportfreq*groupsize, MPI_FLOAT, &status);
		}

	}

	void _consolidate_ordered(int gstep)
	{
		int step = (gstep - reportfreq)/reportfreq;	// normalized step

		numid = 0;	// reset the vector ;-)

//		float all_numbers_unordered[groupsize*reportfreq];	// explicit allocation / vector?
		float *all_numbers_unordered = (float *)malloc(groupsize*reportfreq*sizeof(float));
		MPI_Gather(&numbers[0], reportfreq, MPI_FLOAT, all_numbers_unordered, reportfreq, MPI_FLOAT, 0, group->new_comm);	// 0: master of the group (new_rank == 0)

		int k = 0;
		for (int i = 0; i < reportfreq; i++) {
			for (int j = i; j < groupsize*reportfreq; j+= reportfreq) {
				all_numbers[k] = all_numbers_unordered[j];
				k++;
			}
		}
		free(all_numbers_unordered);

		if (new_rank == 0) {	// master: perform parallel io
			MPI_Status status;

			size_t base = step * (size * reportfreq) * sizeof(float);
			size_t offset  = master_rank* groupsize * sizeof(float);

			MPI_File_set_view(outfile, base+offset, MPI_CHAR, filetype, (char *)"native", MPI_INFO_NULL);
			MPI_File_write_at(outfile, 0, (void *)&all_numbers[0], 1, buftype, &status);

//			int nbytes;
//			MPI_Get_elements(&status, MPI_CHAR, &nbytes);
//			printf("====== number of bytes written = %d ======\n", nbytes); fflush(0);

		}

	}

	void _consolidate_ordered_async(int gstep)
	{
		async_counter++;

/*
		if (lastcall) {
			printf("chunked_async: lastcall = %d\n", lastcall); fflush(0);
		}
*/

		int step = (gstep - reportfreq)/reportfreq;	// normalized step

		if (new_rank == 0) {
			MPI_Status status;
			if (async_counter > 1) MPI_Wait(&request, &status);
		}

		numid = 0;	// reset the vector ;-)

//		float all_numbers_unordered[groupsize*reportfreq];	// explicit allocation / vector?
		float *all_numbers_unordered = (float *)malloc(groupsize*reportfreq*sizeof(float));
		MPI_Gather(&numbers[0], reportfreq, MPI_FLOAT, all_numbers_unordered, reportfreq, MPI_FLOAT, 0, group->new_comm);	// 0: master of the group (new_rank == 0)

		int k = 0;
		for (int i = 0; i < reportfreq; i++) {
			for (int j = i; j < groupsize*reportfreq; j+= reportfreq) {
				all_numbers[k] = all_numbers_unordered[j];
				k++;
			}
		}
		free(all_numbers_unordered);

		if (new_rank == 0) {	// master: perform parallel io
			MPI_Status status;

			size_t base = step * (size * reportfreq) * sizeof(float);
			size_t offset  = master_rank* groupsize * sizeof(float);

			MPI_File_set_view(outfile, base+offset, MPI_CHAR, filetype, (char *)"native", MPI_INFO_NULL);

			MPI_File_iwrite_at(outfile, 0, (void *)&all_numbers[0], 1, buftype, &request);
			if (lastcall) {
				MPI_Wait(&request, &status);
			}


//			int nbytes;
//			MPI_Get_elements(&status, MPI_CHAR, &nbytes);
//			printf("====== number of bytes written = %d ======\n", nbytes); fflush(0);

		}

	}

	void Consolidate(int gstep)
	{
		mygstep = gstep;

		if (layout == 0) {
			_consolidate_chunked_async(gstep);
		}
		else if (layout == 1) {
			_consolidate_ordered(gstep);
		}
		else if (layout == 2) {
			_consolidate_chunked(gstep);
		}
		else if (layout == 3) {
			_consolidate_ordered_async(gstep);
		}
	}

	void Finalize()
	{
		if (numid > 0) {
			lastcall = 1;
			for (int i = numid + 1; i < MAX_REPORTFREQ; i++) numbers[i] = -1;
			Consolidate(mygstep + reportfreq);
		}

		if (new_rank == 0) {	// parallel io
			MPI_File_close(&outfile);
		}

		free(all_numbers);
		free(numbers);
	}

};


#if 0
#include <stdio.h>

using namespace std;

#include "ParIO.h"

int main(int argc, char *argv[])
{
	int rank, size;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	/* Extract the original group handle */
	/* Divide tasks into two distinct groups based upon rank */

	int groupsize = 8;
	int reportfreq = 5;
	int layout = 0;

	MPI_ParIO_Group pio_group;
	MPI_ParIO hist_rank, hist_step;

	pio_group.Init(groupsize, reportfreq, layout, MPI_COMM_WORLD);
	hist_rank.Init((char *)"hist_rank.bin", &pio_group);
	hist_step.Init((char *)"hist_step.bin", &pio_group);

	int step;

	for (step = 0; step < 15; step++) {
		if (step% 5 == 0 && step > 0) {
			hist_rank.Consolidate(step);
			hist_step.Consolidate(step);
		}

//		float number = rank*10.0 + rand()%10  + step*100;
//		float number = rank;
//		float number = step;
		hist_rank.Notify((float)rank);
		hist_step.Notify((float)step);
	}

	hist_rank.Finalize();
	hist_step.Finalize();

	MPI_Finalize();

	return 0;
}
#endif
