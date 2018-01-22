#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "zstd.h"

#include <omp.h>

size_t CompressData( const char* buffIn, size_t read, char* buffOut, size_t buffOutSize )
{
	int cLevel = 1000;
	ZSTD_CStream* const cstream = ZSTD_createCStream();
	size_t const initResult = ZSTD_initCStream(cstream, cLevel);

        ZSTD_inBuffer input = { buffIn, read, 0 };
        ZSTD_outBuffer output = { buffOut, buffOutSize, 0 };
        size_t toRead = ZSTD_compressStream(cstream, &output, &input);

        ZSTD_outBuffer output = { buffOut, buffOutSize, 0 };
	size_t const remainingToFlush = ZSTD_endStream(cstream, &output);   /* close frame */
	if (remainingToFlush) { fprintf(stderr, "not fully flushed"); exit(13); }
	ZSTD_freeCStream(cstream);

	return output.pos;
}

size_t UncompressData( const char* buffIn, size_t buffInSize, char* buffOut, size_t buffOutSize)
{
	ZSTD_DStream* const dstream = ZSTD_createDStream();
	if (dstream==NULL) { fprintf(stderr, "ZSTD_createDStream() error \n"); exit(10); }
	size_t const initResult = ZSTD_initDStream(dstream);
	if (ZSTD_isError(initResult)) { fprintf(stderr, "ZSTD_initDStream() error : %s \n", ZSTD_getErrorName(initResult)); exit(11); }
    
	ZSTD_inBuffer input = { buffIn, buffInSize, 0};
	ZSTD_outBuffer output = { buffOut, buffOutSize, 0};
	size_t toRead = ZSTD_decompressStream(dstream, &output, &input);  /* toRead : size of next compressed block */
	if (ZSTD_isError(toRead)) { fprintf(stderr, "ZSTD_decompressStream() error : %s \n", ZSTD_getErrorName(toRead)); exit(12); }

	ZSTD_freeDStream(dstream);
	return output.pos;
}

#ifndef dtype
#define dtype float
#endif

#ifndef SZ
#define SZ	(32*32*32)
#endif

static dtype indat[SZ];
static dtype outdat[SZ];
static dtype outdat2[SZ];

int main(int argc, char *argv[])
{
//	unsigned char *myout = NULL;
        size_t outs;
        size_t outs2;

        int i;
//	for (i = 0; i < SZ; i++) indat[i] = i;
	srand48(1);
	for (i = 0; i < SZ; i++) indat[i] = drand48();


	double t0 = omp_get_wtime();
	outs = CompressData( (const char*)indat, SZ*sizeof(dtype), (char *) outdat, SZ*sizeof(dtype));
	double t1 = omp_get_wtime();
	
        printf("outs = %ld\n", outs);
        printf("rate = %.2lf\n", (1.0*SZ*sizeof(dtype))/outs);
        printf("time = %.2lf s\n", t1-t0);

	outs2 = UncompressData((const char*)outdat, outs, (char *) outdat2, SZ*sizeof(dtype));
        printf("outs2 = %d\n", outs2);

        for (i = 0; i < SZ; i+=(SZ/4)) printf("%d : %lf vs %lf\n", i, indat[i], outdat2[i]);
        for (i = SZ-3; i < SZ; i+=1) printf("%d : %lf vs %lf\n", i, indat[i], outdat2[i]);

	return 0;
}
