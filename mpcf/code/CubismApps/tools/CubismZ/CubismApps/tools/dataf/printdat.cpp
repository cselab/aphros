#include <stdio.h>
#include <stdlib.h>

#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif

int main(int argc, char *argv[])
{
	long i;
	FILE *fp;
	Real f;

	if (argc != 2) {
		printf("usage: %s <filename>\n", argv[0]);
		exit(1);
	}
	

	fp = fopen(argv[1], "rb");
	if (fp == NULL) {
		printf("fp == NULL\n");
		exit(1);
	}

	while (1)
	{
		int b = fread(&f, 1, sizeof(Real), fp);
		if (b == 0) break;
		
		printf("%.16lf", f);
//		printf("%.4lf ", f);
		i++;
//		if (i % 4 == 0) printf("\n");
		if (i % 1 == 0) printf("\n");


	}

	fclose(fp);	
	return 0;
}
