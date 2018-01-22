#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif

int main(int argc, char *argv[])
{
	int i;
	FILE *fp1, *fp2;
	Real f1;
	double f2;
	if (argc != 3) {
		printf("usage: %s <filename1> <filename2>\n", argv[0]);
		exit(1);
	}
	

	fp1 = fopen(argv[1], "rb");
	if (fp1 == NULL) {
		printf("fp1 == NULL\n");
		exit(1);
	}

	fp2 = fopen(argv[2], "rb");
	if (fp2 == NULL) {
		printf("fp2 == NULL\n");
		exit(1);
	}
	
	while (1)
	{
		int b1 = fread(&f1, 1, sizeof(Real), fp1);
		if (b1 == 0) break;

		int b2 = fread(&f2, 1, sizeof(double), fp2);
		if (b2 == 0) break;
		
		printf("%lf %lf %lf\n", f1, f2, fabs(f1-f2));
	}

	fclose(fp1);	
	fclose(fp2);
	return 0;
}
