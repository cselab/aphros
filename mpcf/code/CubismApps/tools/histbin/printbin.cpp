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

	if (argc != 6) {
		printf("usage: %s <binfile> <nprocs> <nsteps> <repfreq> <nbins>\n", argv[0]);
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
	int Nbins = atoi(argv[5]);


	double minval = +1e9;
	double maxval = -1e9;
	for (int i = 0; i < nprocs*nsteps; i++)
	{
		fread(&f, 1, sizeof(float), fp);
		f = ReverseFloat(f);
		if (f < minval) minval = f;
		if (f > maxval) maxval = f;		
	}
	
	printf("minval = %lf maxval = %lf\n", minval, maxval);

//	double h = (maxval-minval)/nprocs*4;
	double h = (maxval-minval)/16/2;
	double actualmaxval = maxval;
	
	maxval *= 1.1;
	minval = minval - (maxval - actualmaxval);
//	int Nbins = (int) ceil(maxval/h);
//	int Nbins = 10000;

	printf("h = %lf minval = %lf maxval = %lf Nbins = %d\n",  h, minval, maxval, Nbins);

	float meas[nprocs];
	int g_max_nb = 0;

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
		
		double x;
		gsl_histogram *hist = gsl_histogram_alloc (Nbins);
		gsl_histogram_set_ranges_uniform (hist, minval, maxval);

		for (int p = 0; p < nprocs; p++)
			gsl_histogram_increment (hist, meas[p]);

		int max_nb = 0;
		for (int b = 0; b < gsl_histogram_bins(hist); b++) {
			int nb = (int)gsl_histogram_get(hist, b);
			if (nb > max_nb) max_nb = nb;
		}

		gsl_histogram_free (hist);


		if (max_nb > g_max_nb) g_max_nb = max_nb;
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
		
#if 0
		printf("step %3d: ", i);
		for (int p = 0; p < nprocs; p++)
		{
			printf("%.4f ",  meas[p]);
		}
		printf("\n");
#endif

		{
		double x;
		gsl_histogram *hist = gsl_histogram_alloc (Nbins);
		gsl_histogram_set_ranges_uniform (hist, minval, maxval);

		FILE *of;
		char fname[256];

		sprintf(fname, "%s_%04d.txt", argv[1], i);
		of = fopen(fname, "w");
		for (int p = 0; p < nprocs; p++)
			gsl_histogram_increment (hist, meas[p]);

		int max_nb = 0;
		for (int b = 0; b < gsl_histogram_bins(hist); b++) {
			double low, high, med;
			gsl_histogram_get_range(hist, b, &low, &high);
			med = (low + high)/2;
			int nb = (int)gsl_histogram_get(hist, b);
			if (nb > max_nb) max_nb = nb;
			fprintf(of, "%f %d\n",  med, nb); 
		}

		gsl_histogram_free (hist);
		fclose(of);
		
		{
		gnuplot_ctrl * g = NULL;

		g = gnuplot_init();

		gnuplot_cmd(g, "set xrange [%f:%f]", minval, maxval);
//		gnuplot_cmd(g, "set yrange [0:%d]", nprocs);
//		gnuplot_cmd(g, "set yrange [0:100]", nprocs);
		gnuplot_cmd(g, "set yrange [0:%d]", g_max_nb);
#if 0
		gnuplot_cmd(g, "set term x11 %d persist", i);
#else
		gnuplot_cmd(g, "set term png", i);
		gnuplot_cmd(g, "set output \"%s.png\"", fname);
#endif

#if 0
		gnuplot_cmd(g, "plot \"%s\" with lines", fname);
#else
		gnuplot_cmd(g, "set boxwidth %f", h);
		gnuplot_cmd(g, "set style fill solid");
		gnuplot_cmd(g, "plot \"%s\" with boxes", fname);
#endif

		gnuplot_close(g);
		}

		}
	}

	fclose(fp);

	char cmd[256];

	sprintf(cmd, "ffmpeg -y -r 10 -b 1000k -i %s_%%04d.txt.png  %s_movie.avi", argv[1], argv[1]);
	system(cmd);

	return 0;
}
