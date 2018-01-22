// peh
#include <stdio.h>
#include <stdlib.h>

float ReverseFloat( const float inFloat )
{
   float retVal;
   char *floatToConvert = ( char* ) & inFloat;
   char *returnFloat = ( char* ) & retVal;

   // swap the bytes into a temporary buffer
   returnFloat[0] = floatToConvert[3];
   returnFloat[1] = floatToConvert[2];
   returnFloat[2] = floatToConvert[1];
   returnFloat[3] = floatToConvert[0];

   return retVal;
}

int main(int argc, char *argv[])
{
	FILE *fp;
	float f;
	int nprocs = 1;

	fp = fopen(argv[1], "rb");
	if (!fp) {
		printf("file not found!\n");
		exit(1);
	}

	if (argc > 2) {
		nprocs = atoi(argv[2]);
	}

	int i = 0;
	double s = 0;
	while(!feof(fp)) {
		if (fread(&f, 1, sizeof(float), fp) > 0) {
			i++;
#if defined(_SWAP_BYTES_)
			printf("%f ",ReverseFloat(f));
			s += ReverseFloat(f);
#else
			printf("%f ",f);
			s += f;
#endif
			if (i == nprocs) {
				printf("\n");
				i = 0;
			}
		}
	}
	printf("sum = %lf\n", s);
 
	fclose(fp);
}
