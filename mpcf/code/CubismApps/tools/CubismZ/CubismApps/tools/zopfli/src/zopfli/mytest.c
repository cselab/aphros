#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "zopfli.h"
#include <zlib.h>

int zlib_decompress(const void *src, int srcLen, void *dst, int dstLen) {
    z_stream strm  = {0};
    strm.total_in  = strm.avail_in  = srcLen;
    strm.total_out = strm.avail_out = dstLen;
    strm.next_in   = (Bytef *) src;
    strm.next_out  = (Bytef *) dst;

    strm.zalloc = Z_NULL;
    strm.zfree  = Z_NULL;
    strm.opaque = Z_NULL;

    int err = -1;
    int ret = -1;

    err = inflateInit2(&strm, (15 + 32)); //15 window bits, and the +32 tells zlib to to detect if using gzip or zlib
    if (err == Z_OK) {
        err = inflate(&strm, Z_FINISH);
        if (err == Z_STREAM_END) {
            ret = strm.total_out;
        }
        else {
             inflateEnd(&strm);
             return err;
        }
    }
    else {
        inflateEnd(&strm);
        return err;
    }

    inflateEnd(&strm);
    return ret;
}

void shuffle(char *in, int n, int s)
{
	int i, j, k;
	int b;

	char *tmp = malloc(n);

	j = 0;
	for (b = 0; b < s; b++)
	{
		for (i = b; i < n; i+= s) {
			tmp[j] = in[i];
			j++;
		}
	}

	memcpy(in, tmp, n);
	free(tmp);
}

void reshuffle(char *in, int n, int s)
{
	int i, j, k;
	int b;

	char *tmp = malloc(n);
	int c = n/s;
	
	j = 0;
	for (i = 0; i < n/s; i++) {
		for (b = 0; b < s; b++) {
			tmp[j] = in[i + b*c]; j++;

		}
	}

	memcpy(in, tmp, n);
	free(tmp);
}


int main(int argc, char *argv[])
{
#ifndef dtype
#define dtype float
#endif

#define SZ	(32*32*32)
        dtype indat[SZ];
        dtype outdat[SZ];
	unsigned char *myout = NULL;
        size_t outs;
        dtype outdat2[SZ];
        size_t outs2;

        int i;
//	for (i = 0; i < SZ; i++) indat[i] = i;
	srand48(1);
	for (i = 0; i < SZ; i++) indat[i] = drand48();

	shuffle((char *)indat, SZ*sizeof(dtype), sizeof(dtype));

	ZopfliOptions options;
	ZopfliInitOptions(&options);


	ZopfliCompress(&options, ZOPFLI_FORMAT_ZLIB,
		(const unsigned char *) indat, SZ*sizeof(dtype),
		(unsigned char **) &myout, (size_t *) &outs);

	
        printf("outs = %ld\n", outs);
        printf("rate = %.2lf\n", (1.0*SZ*sizeof(dtype))/outs);

	outs2 = zlib_decompress((const void *)myout, outs, (void *) outdat2, SZ*sizeof(dtype));
//	outs2 = SZ*sizeof(dtype);

//	outs2 = LZ4_uncompress_unknownOutputSize ((const char*) outdat, (char *)outdat2, outs, SZ*sizeof(dtype));
        printf("outs2 = %d\n", outs2);
	reshuffle((char *)outdat2, SZ*sizeof(dtype), sizeof(dtype));

        for (i = 0; i < SZ; i+=1024) printf("%d : %lf\n", i, outdat2[i]);
        for (i = SZ-10; i < SZ; i+=1) printf("%d : %lf\n", i, outdat2[i]);

 return 0;
}

