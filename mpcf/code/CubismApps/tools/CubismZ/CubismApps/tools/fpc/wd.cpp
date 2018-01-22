#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif

#define _BLOCKSIZE_	32

int main(int argc, char *argv[])
{
	int i;
	FILE *fp;
	Real f;


	if (argc != 4) {
		printf("usage: %s <filename> <n> <seed>\n", argv[0]);
		exit(1);
	}

	int n = atoi(argv[2]);
	int seed = atoi(argv[3]);
	
	srand48(seed);

	fp = fopen(argv[1], "wb");
	if (fp == NULL) {
		printf("fp == NULL\n");
		exit(1);
	}
	
	int s = n; 
	for (i = 0; i < s; i++) {
		f = i;
		fwrite(&f, 1, sizeof(Real), fp);
	}

	printf("datasize = %ld  Bytes\n", s*sizeof(Real));
	printf("datasize = %ld KBbytes\n", s*sizeof(Real)/(1024));
	printf("datasize = %ld MBytes\n", s*sizeof(Real)/(1024*1024));
	fclose(fp);	
	return 0;
}
