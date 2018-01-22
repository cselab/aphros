#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{
	FILE *fp1, *fp2;
	FILE *fp3;

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

	char name3[256];
	sprintf(name3, "%s_diff.bin", argv[1]);
	fp3 = fopen(name3, "wb");

	long n = 0;
	float e_inf = 0;
	float e_1 = 0;
	float e_2 = 0;
	
	float n_inf = 0;
	float n_1 = 0;
	float n_2 = 0;

	float threshold = 0.0001; //0.1;	// 0.1%

	//printf("give threshold =");
	//scanf("%f", &threshold);

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
			n++;

#if 1	// absolute error
			if (err >= threshold)
			{
				over_counter++;
				//printf("%5d: %15.8f %15.8f (%15.8f)\n", over_counter, f1, f2, err); 
				err += 1.0;
				fwrite(&err, 1, 4, fp3);
				//fwrite(&one, 1, 4, fp3);
			}
			else
			{
				float zero = 0.0;
				fwrite(&zero, 1, 4, fp3);
			}
#else
			float rel_err = 100.0 * (err / v);
			if (rel_err >= threshold)
			{
				//float one = 1;
				over_counter++;
				//printf("%5d: %15.8f %15.8f (%15.8f)\n", over_counter, f1, f2, rel_err); 
				err += 1.0;
				fwrite(&err, 1, 4, fp3);
				//fwrite(&one, 1, 4, fp3);
			}
			else
			{
				float zero = 0.0;
				fwrite(&zero, 1, 4, fp3);
			}
#endif

		}

		if ((r1 == 0)&&(r2 == 0))
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
	printf("rel(e_1)   = %.8f\n", e_1/n_1);
//	printf("mean(e_1)  = %f\n", e_1/n);
//	printf("\n");

	printf("e_2        = %.8f\n", sqrt(e_2));
	printf("rel(e_2)   = %.8f\n", sqrt(e_2)/sqrt(n_2));
//	printf("mean(e_2)  = %f\n", sqrt(e_2)/n);
	printf("\n");

	printf("counter    = %ld\n", over_counter);
	printf("percentage = %.8f\n", (float)over_counter/(float)n);

	fclose(fp1);
	fclose(fp2);
	fclose(fp3);

	return 0;
}
