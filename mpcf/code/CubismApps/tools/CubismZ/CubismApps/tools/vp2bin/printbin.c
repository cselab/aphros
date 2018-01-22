#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{
	FILE *fp1;

	if (argc != 2) {
		printf("usage: %s <filename1>\n", argv[0]);
		exit(1);
	}

	fp1 = fopen(argv[1], "r");
	if ((fp1 == NULL)) {
		printf("file not found!\n");
		exit(1);
	}

	long n = 0;
	while (1)
	{
		float f1;
		int r1;

		r1 = fread(&f1, 1, 4, fp1);

		if ((r1 == 4))
		{
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

			printf("%15.8f @ [%4d,%5d]:(%2d,%2d,%2d)\n", f1, blockid, elemid, x, y, z); 
#endif
			n++;
		}

		if ((r1 == 0))
		{
			printf("finishing gracefully\n");
			break;
		}

	}

	printf("=========================\n");
	fclose(fp1);

	return 0;
}
