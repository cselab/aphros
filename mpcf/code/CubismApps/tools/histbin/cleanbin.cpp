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
	FILE *fp, *ofp;
	float f;
	int nprocs = 1;

	fp = fopen(argv[1], "rb");
	if (!fp) {
		printf("file not found!\n");
		exit(1);
	}

	ofp = fopen(argv[2], "wb");
	if (!ofp) {
		printf("file was not created!\n");
		exit(1);
	}

	int i = 0;
	while(!feof(fp)) {
		if (fread(&f, 1, sizeof(float), fp) > 0) {
			i++;
#if defined(_SWAP_BYTES_)
			printf("%f ",ReverseFloat(f));
#else
			printf("%f ",f);
#endif
			if (i == nprocs) {
				printf("\n");
				i = 0;
			}

#if defined(_SWAP_BYTES_)
			if (ReverseFloat(f) >= 0.0)
#else
			if (f >= 0.0)
#endif
			{
				fwrite(&f, 1, sizeof(float), ofp);
			}
			
		}
	}
 
	fclose(fp);
	fclose(ofp);
	
	return 0;
}
