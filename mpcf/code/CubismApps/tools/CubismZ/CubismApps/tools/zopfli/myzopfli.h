#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "zopfli.h"
#include <zlib.h>

static int zopfli_compress(const unsigned char * indat, size_t inbytes, unsigned char * outdat)
{
	unsigned char *myout = 0; 
	size_t outs = 0;

	ZopfliOptions options;
	ZopfliInitOptions(&options);
//	options.numiterations = 1;
	ZopfliCompress(&options, ZOPFLI_FORMAT_ZLIB,
		(const unsigned char *) indat, inbytes,
		(unsigned char **) &myout, (size_t *) &outs);
	memcpy(outdat, myout, outs);
	free(myout); 
	return outs;
}

static int zopfli_decompress(const void *src, int srcLen, void *dst, int dstLen) {
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
