// peh
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <gsl/gsl_histogram.h>
extern "C"
{
#include "gnuplot_i.h"
}

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

	if (argc != 2) {
		printf("usage: %s <nsteps>\n", argv[0]);
		exit(1);
	}

	int nsteps = atoi(argv[1]);

#if 0
	double max_val = -1e8;
	for (int i = 0; i < nsteps; i++)
	{
		double val;
		char fname[256];
		sprintf(fname, "lab_load.%04d.txt", i);
		
		FILE *fp = fopen(fname, "r");
		if (fp == NULL) {
			printf("file <%s> not found!\n", fname);
			exit(1);
		}
		while (fscanf(fp, "%lf", &val) != EOF) {
			if (val > max_val) max_val = val;	
		}
		fclose(fp);
	}
#else
	double max_val = 2*1e-3; // 2 ms
#endif


	for (int i = 0; i < nsteps; i++)
	{
		char fname[256];
		sprintf(fname, "lab_load.%04d.txt", i);
		
		FILE *fp = fopen(fname, "r");
		if (fp == NULL) {
			printf("file <%s> not found!\n", fname);
			exit(1);
		}
		fclose(fp);
		
		{
		gnuplot_ctrl * g = NULL;

		g = gnuplot_init();

//		gnuplot_cmd(g, "set xrange [%f:%f]", minval, maxval);
//		gnuplot_cmd(g, "set yrange [0:%d]", nprocs);
//		gnuplot_cmd(g, "set yrange [0:100]", nprocs);
//		gnuplot_cmd(g, "set yrange [0:%d]", g_max_nb);
		gnuplot_cmd(g, "set yrange [0:%f]", max_val);
#if 0
		gnuplot_cmd(g, "set term x11 %d persist", i);
#else
		gnuplot_cmd(g, "set term png", i);
		gnuplot_cmd(g, "set output \"%s.png\"", fname);
#endif

#if 1
		gnuplot_cmd(g, "plot \"%s\" with lines", fname);
#else
		gnuplot_cmd(g, "set boxwidth %f", h);
		gnuplot_cmd(g, "set style fill solid");
		gnuplot_cmd(g, "plot \"%s\" with boxes", fname);
#endif

		gnuplot_close(g);
		}

	}


	char cmd[256];

	sprintf(cmd, "ffmpeg -y -r 10 -b 1000k -i lab_load.%%04d.txt.png  lab_load_movie.avi");
	system(cmd);

	return 0;
}
