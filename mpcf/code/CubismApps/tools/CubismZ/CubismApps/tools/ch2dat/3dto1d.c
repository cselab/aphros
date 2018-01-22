#include <stdio.h>


int main(int argc, char *argv[])
{
	int bs = 64;

	int zLength, yLength;
	zLength = yLength = bs;
	int i;

//	i = atoi(argv[1]);

	while (scanf("%d", &i) > 0)
	{
		int iz = i % zLength;
		int iy = (i / zLength) % yLength;
		int ix = i / (yLength * zLength); 

		printf("%d -> (%d,%d,%d)\n", i, ix, iy, iz);
	}
	
	return 0;
}
