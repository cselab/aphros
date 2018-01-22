#include <stdio.h>
#include <math.h>

int gsize = 32;

void pos(float p[3], int ix, int iy, int iz)
{
	float gextend = 1.0*gsize;

	p[0] = ix/gextend;
	p[1] = iy/gextend;
	p[2] = iz/gextend;
}

int main(int argc, char *argv[])
{
	int bsize = 32;

	FILE *fp = fopen("in.dat", "wb");

	for(int iz=0; iz<bsize; iz++)
	for(int iy=0; iy<bsize; iy++)
	for(int ix=0; ix<bsize; ix++)
	{
		float p[3];
		pos(p, ix, iy, iz);
		float d = sqrt(pow(p[0]-0.5,2)+pow(p[1]-0.5,2)+pow(p[2]-0.5,2))-0.1;
		fwrite (&d, 1, sizeof(float),fp);
		//printf("%d %d %d -> %f\n", ix, iy, iz, d);

	}
	fclose(fp);


	return 0;
}

