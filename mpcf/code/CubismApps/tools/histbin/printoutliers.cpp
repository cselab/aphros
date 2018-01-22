// peh
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <gsl/gsl_histogram.h>

#include "gnuplot_i.h"

float ReverseFloat( const float inFloat )
{
	float retVal;
#if defined(_SWAP_BYTES_)
	char *floatToConvert = ( char* ) & inFloat;
	char *returnFloat = ( char* ) & retVal;

	// swap the bytes into a temporary buffer
	returnFloat[0] = floatToConvert[3];
	returnFloat[1] = floatToConvert[2];
	returnFloat[2] = floatToConvert[1];
	returnFloat[3] = floatToConvert[0];
#else
	retVal = inFloat;
#endif
	return retVal;
}

int main(int argc, char *argv[])
{
	FILE *fp;
	float f;

	if (argc != 5) {
		printf("usage: %s <binfile> <nprocs> <nsteps> <repfreq>\n", argv[0]);
		exit(1);
	}

	fp = fopen(argv[1], "rb");
	if (!fp) {
		printf("file not found!\n");
		exit(1);
	}

	int nprocs = atoi(argv[2]);
	int nsteps = atoi(argv[3]);
	int repfreq = atoi(argv[4]);


	double minval = +1e9;
	double maxval = -1e9;
	double sumval = 0.0;
	double avgval;
	for (int i = 0; i < nprocs*nsteps; i++)
	{
		fread(&f, 1, sizeof(float), fp);
		f = ReverseFloat(f);
		sumval += f;
		
		if (f < minval) minval = f;
		if (f > maxval) maxval = f;
	}
	
	avgval = sumval / (nprocs*nsteps);

	fseek(fp, 0, SEEK_SET);
	int N = nprocs*nsteps;
	double sumstd = 0.0;
	for (int i = 0; i < N; i++)
	{
		fread(&f, 1, sizeof(float), fp);
		f = ReverseFloat(f);
		sumstd += pow(f - avgval, 2);
	}

	double stdval = sqrt(sumstd/(N-1));

	printf("minval  = %lf\n", minval);
	printf("maxval  = %lf\n", maxval);
	printf("avgval  = %lf\n", avgval);
	printf("stdval1 = %lf\n", stdval);
	printf("std/avg = %.2lf%%\n", 100.0*stdval/avgval);

	float meas[nprocs];
	int g_max_nb = 0;


	int **outliers;
	outliers = (int **)malloc(nsteps*sizeof(int *));
	for (int i = 0; i < nsteps; i++) {
		outliers[i] = (int *)malloc(nprocs*sizeof(int));
		for (int p = 0; p < nprocs; p++)
			outliers[i][p] = -1;
	}

	for (int i = 0; i < nsteps; i++)
	{
	        int posvec[nprocs];
                int base = (i / repfreq) * (repfreq*nprocs) + i % repfreq;
                for (int p = 0; p < nprocs; p++)
                {
                        int pos = base + p*repfreq;
                        posvec[p] = pos;
                }

                for (int p = 0; p < nprocs; p++)
		{
			fseek(fp, posvec[p]*sizeof(float), SEEK_SET);
			fread(&f, 1, sizeof(float), fp);
			meas[p] = ReverseFloat(f);
		}

		int pos = 0;
                for (int p = 0; p < nprocs; p++)
		{
//			if (fabs(meas[p]-avgval) > 5*stdval)
			if ((meas[p]-avgval) > 3*stdval)
			{
				printf("step=%4d  proc=%4d time=%lf\n", i, p, meas[p]); 
				
				outliers[i][pos] = p; pos++;
			}
		}
		
	}


	for (int i = 0; i < nsteps; i++)
	{
		printf("step %4d:", i);
		for (int pos = 0; pos < nprocs; pos++)
		{
			if (outliers[i][pos] == -1) break;
			printf("%4d ", outliers[i][pos]);
		}
		printf("\n");
	}


	int times_per_proc[nprocs];
	for (int p = 0; p < nprocs; p++) times_per_proc[p] = 0;
	
	for (int i = 0; i < nsteps; i++)
	{
		for (int pos = 0; pos < nprocs; pos++)
		{
			if (outliers[i][pos] == -1) break;
			times_per_proc[outliers[i][pos]]++;
		}
	}
	
	printf("TIMES PER PROC\n");
	for (int p = 0; p < nprocs; p++)
	{
		//if (times_per_proc[p] != 0)
		{
			printf("%4d : %4d\n", p, times_per_proc[p]);
		}
	}	


	fclose(fp);
}
