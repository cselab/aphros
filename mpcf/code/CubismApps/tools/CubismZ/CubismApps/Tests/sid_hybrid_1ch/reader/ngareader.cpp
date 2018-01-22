/*
 *  ngareader.cpp
 *
 *  Created by Panos Hadjidoukas on 4/4/16.
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
	int full; // = 1;
	int z_idx; // = 0;
	int y_idx; // = 0;

	if (argc != 4) {
		printf("usage: %s <z> <y> <full:0 or 1>\n", argv[0]);
		exit(1);
	}

	z_idx = atoi(argv[1]);
	y_idx = atoi(argv[2]);
	full  = atoi(argv[3]);



	char *filename = (char *)"../data/data.init";

	///////////////////////////////////////

	FILE *fid = fopen(filename,"rb");
	int dims[4];

	fread(dims, sizeof(int), 4, fid);	// read only the first 4 elements from the file

	// Dimensions 
	int nx      = dims[0];
	int ny      = dims[1];
	int nz      = dims[2];
	int nvar    = dims[3];


	printf("nx ny nz = %d %d %d\n", nx, ny, nz);

	int nsize   = nx*ny*nz;

	//  The following lines are necessary to move the file pointer to the correct
	// location befory you start reading in the variables

	double dt, time;

	fread(&dt, sizeof(double), 1, fid);	//skipping over 'dt'
	fread(&time, sizeof(double), 1, fid);	//skipping over 'time'

	char varnames[8][nvar];

	for (int var=0; var<nvar; var++) {
		fread(varnames[var], sizeof(char), 8, fid);	//skipping over the stored names
	}

	for (int var=0; var<nvar; var++) {
		printf("var[%d] = ->%s<-\n", var, varnames[var]);
	}

	// There are usually 5 variables (or nvar number of variables) stored in the data
	// file: U, V, W, P, ZMIX. Each fread call makes the file pointer shift to the end 
	// of the number of values read. Hence each subsequent fread call will start after 
	// where the previous fread call stopped. Add on extra variables if necessary.

	double *dummy;
	double *U, *V, *W, *P;

	//dummy = malloc(nsize*sizeof(double));
	U = (double *)malloc(nsize*sizeof(double));
	V = (double *)malloc(nsize*sizeof(double));
	W = (double *)malloc(nsize*sizeof(double));
	P = (double *)malloc(nsize*sizeof(double));


	// dummy   = fread(fid,nsize,'real*8','ieee-le');    % will produce a column vector
	//U       = reshape(dummy,nx,ny,nz);                % now turning the column vector into a 3D matrix

	//fread(dummy, sizeof(double), nsize, fid);
	if (full == 1) {
		fread(U, sizeof(double), nsize, fid);
	}
	else {
		int z = z_idx;
		int y = y_idx;
		int varid = 0;
		long header = 4*sizeof(int) + 2*sizeof(double) + 8*nvar*sizeof(char); 
		long base = nx*nx*ny*sizeof(double)*varid;
		long offset = nx*ny*z*sizeof(double) + nx*y*sizeof(double);
		fseek(fid, header + base + offset, SEEK_SET);
		//fread(U, sizeof(double), nx*ny, fid);
		fread(U, sizeof(double), nx, fid);
		 
		//fread(&U[0], sizeof(double), nsize/2, fid);
		//fread(&U[nsize/2], sizeof(double), nsize, fid);
	}

	fclose(fid);

	// Now use U, V, W, P, ZMIX for computations you need

	if (full == 1) {
		for (int z = z_idx; z <= z_idx; z++) {
			//for (int y = 0; y < ny; y+=16)
			for (int y = y_idx; y <= y_idx; y++) {
				for (int x = 0; x < nx; x+=16)
					printf("Ux0[%03d,%03d,%03d]=U[%06d]= %f\n", x, y, z, z*nx*ny+y*nx+x, U[z*nx*ny + y*nx + x]); 
			}
		}
	}
	else {
		for (int z = z_idx; z <= z_idx; z++) {
			int l_z = z_idx - z_idx;
			//for (int y = 0; y < ny; y+=16) {
			for (int y = y_idx; y <= y_idx; y++) {
				int l_y = y - y_idx;
				for (int x = 0; x < nx; x+=16)
					//printf("Ux1[%03d,%03d,%03d]=U[%06d]= %f\n", x, y, z, z*nx*ny+y*nx+x, U[l_z*nx*ny + y*nx + x]); 
					printf("Ux1[%03d,%03d,%03d]=U[%06d]= %f\n", x, y, z, z*nx*ny+y*nx+x, U[l_z*nx*ny + l_y*nx + x]); 
			}
		}
	}


	return 0;
}

