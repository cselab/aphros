#include <stdio.h>
#include <stdlib.h>

#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif

int main(int argc, char *argv[])
{
	int i;
	FILE *fp;
	Real f;
	Real dummy;
	int stride = 0;
	int skip = 0;

	if (argc != 4) {
		printf("usage: %s <filename> <skip> <stride>\n", argv[0]);
		exit(1);
	}
	

	fp = fopen(argv[1], "rb");
	if (fp == NULL) {
		printf("fp == NULL\n");
		exit(1);
	}

	skip = atoi(argv[2]);
	for (i = 0; i < skip; i++) {
		int b = fread(&dummy, 1, sizeof(Real), fp);
	}

	stride = atoi(argv[3])-1;
	
	while (1)
	{
		int b = fread(&f, 1, sizeof(Real), fp);
		if (b == 0) break;
		for (i = 0; i < stride; i++)
			b = fread(&dummy, 1, sizeof(Real), fp);
		
		printf("%.16f\n", f);
	}

	fclose(fp);	
	return 0;
}
