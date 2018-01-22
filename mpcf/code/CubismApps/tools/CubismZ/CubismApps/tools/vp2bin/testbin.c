#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{
	FILE *fp1, *fp2;

	if (argc != 3) {
		printf("usage: %s <filename1> <filename2>\n", argv[0]);
		exit(1);
	}

	fp1 = fopen(argv[1], "r");
	fp2 = fopen(argv[2], "r");
	if ((fp1 == NULL)||(fp2 == NULL)) {
		printf("at least one file not found!\n");
		exit(1);
	}

	long n = 0;
	float e_inf = 0;
	float e_1 = 0;
	float e_2 = 0;
	
	float n_inf = 0;
	float n_1 = 0;
	float n_2 = 0;

	float threshold = 0.1;	// 0.1%

	printf("give threshold =");
	scanf("%f", &threshold);

	long over_counter = 0;

	while (1)
	{
		float f1, f2;
		int r1, r2;

		r1 = fread(&f1, 1, 4, fp1);
		r2 = fread(&f2, 1, 4, fp2);

		if (r1 != r2) 
		{
			printf("this should not happen!\n");
			exit(1);
		}

		if ((r1 == 4)&&(r2 == 4))
		{
			float v = fabs(f2);
			if (v > n_inf) n_inf = v;
			n_1 += v;
			n_2 += v*v;

			float err = fabs(f1-f2);

			if (err > e_inf) e_inf = err;
			e_1 += err;
			e_2 += err*err;
			n++;	// number id

			float rel_err;

			if (v != 0.0)
				rel_err = 100.0 * (err / v);
			else
				rel_err = 0.0;

#if 0	// absolute error
			if (err >= threshold)
#else	// absolute relative % error
			if (rel_err >= threshold)
//			if ((rel_err >= threshold)&&(err >= 0.001)) //&&(v > 0.001))
#endif
			{
				over_counter++;
				//printf("%5d: %15.8f %15.8f (%15.8f)\n", over_counter, f1, f2, err); 
				//printf("%5d: %15.8f %15.8f (%15.8f)\n", over_counter, f1, f2, rel_err); 

#if 1
				#ifndef BS
				#define BS	(32)
				#endif
				#define BS2	(BS*BS)
				#define BS3	(BS*BS*BS)
				int blockid = n / BS3;
				int elemid = n % (BS3);

				int x = elemid / (BS2);
				int y = (elemid / BS) % BS;
				int z = elemid % BS;

				printf("%5d: %15.8f %15.8f (%15.8f - %10.2f%%) [%4d,%5d]:(%2d,%2d,%2d)\n", over_counter, f1, f2, err, rel_err, blockid, elemid, x, y, z); 
#endif
			}
		}

		if ((r1 == 0)&&(r2 == 0))
		{
			printf("finishing gracefully\n");
			break;
		}

	}

	printf("=========================\n");
	printf("e_inf      = %.6f\n", e_inf);
//	printf("rel(e_inf) = %f\n", e_inf/n_inf);
//	printf("\n");

	printf("e_1        = %.6f\n", e_1);
	printf("rel(e_1)   = %.3f\n", e_1/n_1);
//	printf("mean(e_1)  = %f\n", e_1/n);
//	printf("\n");

	printf("e_2        = %.6f\n", sqrt(e_2));
	printf("rel(e_2)   = %.6f\n", sqrt(e_2)/sqrt(n_2));
//	printf("mean(e_2)  = %f\n", sqrt(e_2)/n);
	printf("\n");

	printf("counter    = %ld\n", over_counter);
	printf("percentage = %.3f%%\n", 100.0*(float)over_counter/(float)n);

	fclose(fp1);
	fclose(fp2);

	return 0;
}
