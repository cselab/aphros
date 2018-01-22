#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif

int main(int argc, char *argv[])
{
	FILE *fp1, *fp2;

	if ((argc != 3)&&(argc != 4)) {
		printf("usage: %s <filename1> <filename2> [<compressedfilename1>]\n", argv[0]);
		exit(1);
	}

	fp1 = fopen(argv[1], "r");
	fp2 = fopen(argv[2], "r");
	if ((fp1 == NULL)||(fp2 == NULL)) {
		printf("at least one file not found!\n");
		exit(1);
	}


	long n = 0;
	double e_inf = 0;
	double e_1 = 0;
	long double e_2 = 0;
	
	double n_inf = 0;
	double n_1 = 0;
	long double n_2 = 0;

	double maxdata = -DBL_MAX;
	double mindata =  DBL_MAX;

	while (1)
	{
		Real f1, f2;	// f2: reference
		int r1, r2;

		r1 = fread(&f1, 1, sizeof(Real), fp1);
		r2 = fread(&f2, 1, sizeof(Real), fp2);

		if (r1 != r2) 
		{
			printf("this should not happen!\n");
			exit(1);
		}

		if ((r1 == sizeof(Real))&&(r2 == sizeof(Real)))
		{
			Real v = fabs(f2);

			if ((double)v > maxdata) maxdata = (double)v;
			if ((double)v < mindata) mindata = (double)v;

			if (v > n_inf) {
				n_inf = v;
			}
			n_1 += v;
			n_2 += v*v;

			const double err = fabs((double)f1-(double)f2);

			if (err > e_inf) {
				e_inf = err;
				printf("%15.8f vs %15.8f -> %15.8f (rel_err = %15.8f%%)\n", f1, f2, err, 100.0*(err/v));
			}

			e_1 += err;
			e_2 += err*err;
			n++;
		}

		if ((r1 == 0)&&(r2 == 0))
		{
			printf("finishing gracefully\n");
			break;
		}

	}

        fclose(fp1);
        fclose(fp2);


	printf("=========================\n");
	printf("nall 	   = %ld\n", n);
        printf("e_inf	   = %.8f\n", e_inf);
//	printf("rel(e_inf) = %f\n", e_inf/n_inf);
//	printf("\n");

        printf("e_1        = %.8f\n", e_1);
        printf("rel(e_1)   = %.8f\n", e_1/n_1);
	printf("mean(e_1)  = %.6f\n", e_1/n);
//	printf("\n");

        printf("e_2        = %.8f\n", sqrt(e_2));
        printf("rel(e_2)   = %.8f\n", sqrt(e_2)/sqrt(n_2));
	printf("mean(e_2)  = %.6f\n", sqrt(e_2)/n);
        printf("\n");


	double linf = e_inf;
	long nall = n;
	double l1 = e_1 / nall;
	const long double mse = e_2 / nall;
	double l2 = sqrt(e_2) / nall;


	double uncompressed_footprint = sizeof(Real) * nall; 
	double compressed_footprint = uncompressed_footprint; 

#if 1
	if (argc == 4)
	{
		FILE *fpz = fopen(argv[3], "rb");
		if (fpz == NULL)
		{
			printf("fp == NULL for %s\n", argv[3]);
                }
		else
		{
			printf("opening %s\n", argv[3]);
			int fd = fileno(fpz); //if you have a stream (e.g. from fopen), not a file descriptor.
			struct stat buf;
			fstat(fd, &buf);
			compressed_footprint = buf.st_size;
			fclose(fpz);
		}
	}
#endif


	printf("compression-rate: %.2f rel-linf-error: %e rel-mean-error: %e\n",
           uncompressed_footprint / compressed_footprint,
           linf, l1);

	//https://en.wikipedia.org/wiki/Peak_signal-to-noise_ratio
	printf("mindata = %f\n", mindata);
	printf("maxdata = %f\n", maxdata);

	double psnr;

//	if (normalize)
//		psnr = 10 * (double)log10(1. * 1. / mse);
//	else
		psnr = 20 * log10((maxdata - mindata) / (2 * sqrt(mse)));


	printf("PSNR-vs-BITRATE: %.04f bps %.04f dB\n",
		compressed_footprint * 8. / nall, psnr);


}



#if 0
void performance()
{
    printf("assessing performance...\n");

    double linf = 0, l1 = 0;
    long double l2 = 0;

    for(int i = 0; i < nall; ++i)
    {
     	const double e = fabs((double)reference_data[i] - (double)result_data[i]);

        linf = fmax(linf, e);
        l1 += e;
        l2 += e * e;
    }

    l1 /= nall;
    const long double mse = l2 / nall;
    l2 = sqrt(l2) / nall;

    printf("compression-rate: %.2f rel-linf-error: %e rel-mean-error: %e\n",
           sizeof(float) * nall * 1. / compressed_footprint,
           linf, l1);

    //https://en.wikipedia.org/wiki/Peak_signal-to-noise_ratio
    double psnr;

    if (normalize)
        psnr = 10 * (double)log10(1. * 1. / mse);
    else
        psnr = 20 * log10((maxdata - mindata) / (2 * sqrt(mse)));


    printf("PSNR-vs-BITRATE: %.04f bps %.04f dB\n",
           compressed_footprint * 8. / nall, psnr);	// 8x from bytes to bits
}
#endif
