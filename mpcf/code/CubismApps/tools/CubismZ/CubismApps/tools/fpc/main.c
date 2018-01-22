#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "fpc.h"

int main(int argc, char *argv[])
{
 long val, ioc;

 assert(4 <= sizeof(long));
 assert(8 == sizeof(long long));
// assert(0 < SIZE);
// assert(0 == (SIZE & 0xf));
 val = 1;
 assert(1 == *((char *)&val));

#define SZ	(32*32*32)
	double indat[SZ];
	double outdat[SZ];
	int outs;
	double outdat2[SZ];
	int outs2;
	
	int i;
	for (i = 0; i < SZ; i++) indat[i] = i;
	fpc_compress((char *)indat, SZ*8, (char *)outdat, &outs, 10);
	printf("outs = %d\n", outs);
	printf("rate = %.2lf\n", (SZ*8.0)/outs);

	fpc_decompress((char *)outdat, outs, (char *)outdat2, &outs2);
	printf("outs2 = %d\n", outs2);
	
	for (i = 0; i < SZ; i+=1024) printf("%d : %lf\n", i, outdat2[i]);
	for (i = SZ-10; i < SZ; i+=1) printf("%d : %lf\n", i, outdat2[i]);	

 return 0;
}

