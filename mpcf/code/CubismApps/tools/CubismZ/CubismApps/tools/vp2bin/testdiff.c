#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{
	FILE *fp1, *fp2;

	if (argc != 2) {
		printf("usage: %s <filename1>\n", argv[0]);
		exit(1);
	}

	fp1 = fopen(argv[1], "r");
	if ((fp1 == NULL)) {
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
		float f1;
		int r1;

		r1 = fread(&f1, 1, 4, fp1);

		if ((r1 == 4))
		{
			float err = fabs(f1);

			if (err > e_inf) e_inf = err;
			e_1 += err;
			e_2 += err*err;
			n++;

#if 1	// absolute error
			if (err >= threshold)
			{
				over_counter++;
				printf("%5d: %15.8f (%15.8f)\n", over_counter, f1, err); 
			}
#else
			float rel_err = 100.0 * (err);
			if (rel_err >= threshold)
			{
				over_counter++;
				//printf("%5d: %15.8f %15.8f (%15.8f)\n", over_counter, f1, f2, rel_err); 
			}
#endif

		}

		if (r1 == 0)
		{
			printf("finishing gracefully\n");
			break;
		}

	}

	printf("=========================\n");
	printf("e_inf      = %.8f\n", e_inf);
//	printf("rel(e_inf) = %f\n", e_inf/n_inf);
//	printf("\n");

	printf("e_1        = %.8f\n", e_1);
//	printf("rel(e_1)   = %.8f\n", e_1/n_1);
//	printf("mean(e_1)  = %f\n", e_1/n);
//	printf("\n");

	printf("e_2        = %.8f\n", sqrt(e_2));
//	printf("rel(e_2)   = %.8f\n", sqrt(e_2)/sqrt(n_2));
//	printf("mean(e_2)  = %f\n", sqrt(e_2)/n);
	printf("\n");

	printf("counter    = %ld\n", over_counter);
	printf("percentage = %.8f\n", (float)over_counter/(float)n);

	fclose(fp1);

	return 0;
}
