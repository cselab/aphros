/**
 *  @file test_zlib.c
 *  @author Sheng Di
 *  @date June, 2016
 *  @brief gzip compressor code: the interface to call zlib
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <zutil.h>
#include <sz.h>

#define CHECK_ERR(err, msg) { \
    if (err != Z_OK && err != Z_STREAM_END) { \
        fprintf(stderr, "%s error: %d\n", msg, err); \
        exit(1); \
    } \
}

/*zlib_compress() is only valid for median-size data compression. */
unsigned long zlib_compress(unsigned char* data, unsigned long dataLength, unsigned char** compressBytes, int level)
{
	unsigned long outSize;
	unsigned long* outSize_ = (unsigned long*)malloc(sizeof(unsigned long));
	*outSize_ = 0;
	*outSize_ = dataLength;
	*compressBytes = (unsigned char*)malloc(sizeof(unsigned char)*dataLength);
	int err = compress2(*compressBytes, outSize_, data, dataLength, level);
	printf("err=%d\n", err);
	outSize = *outSize_;
	free(outSize_);
	return outSize;
}

unsigned long zlib_compress2(unsigned char* data, unsigned long dataLength, unsigned char** compressBytes, int level)
{
	unsigned long outSize;
	*compressBytes = (unsigned char*)malloc(sizeof(unsigned char)*dataLength);
	
	z_stream stream;
    int err;

    stream.next_in = data;
    stream.avail_in = dataLength;
#ifdef MAXSEG_64K
    /* Check for source > 64K on 16-bit machine: */
    if ((uLong)stream.avail_in != dataLength) return Z_BUF_ERROR;
#endif
    stream.next_out = *compressBytes;
    stream.avail_out = dataLength;
    if ((uLong)stream.avail_out != dataLength) return Z_BUF_ERROR;

    stream.zalloc = (alloc_func)0;
    stream.zfree = (free_func)0;
    stream.opaque = (voidpf)0;
//	stream.data_type = Z_TEXT;

    //err = deflateInit(&stream, level); //default  windowBits == 15.
    int windowBits = 14; //8-15
    if(szMode==SZ_BEST_COMPRESSION)
		windowBits = 15;
	
    err = deflateInit2(&stream, level, Z_DEFLATED, windowBits, DEF_MEM_LEVEL,
                         Z_DEFAULT_STRATEGY);//Z_FIXED); //Z_DEFAULT_STRATEGY
    if (err != Z_OK) return err;

    err = deflate(&stream, Z_FINISH);
    if (err != Z_STREAM_END) {
        deflateEnd(&stream);
        return err == Z_OK ? Z_BUF_ERROR : err;
    }

    err = deflateEnd(&stream);
    
    outSize = stream.total_out;
    return outSize;
}

unsigned long zlib_uncompress(unsigned char* compressBytes, unsigned long cmpSize, unsigned char** oriData, unsigned long targetOriSize)
{
	unsigned long outSize;
	unsigned long* outSize_ = (unsigned long*)malloc(sizeof(unsigned long));
	*oriData = (unsigned char*)malloc(sizeof(unsigned char)*targetOriSize);
	uncompress(*oriData, outSize_, compressBytes, cmpSize); 
	outSize = *outSize_;
	free(outSize_);
	return outSize;
}

unsigned long zlib_uncompress2 (unsigned char* compressBytes, unsigned long cmpSize, unsigned char** oriData, unsigned long targetOriSize)
{
    z_stream stream;

	unsigned long outSize;
	*oriData = (unsigned char*)malloc(sizeof(unsigned char)*targetOriSize);

    stream.zalloc = Z_NULL;
    stream.zfree = Z_NULL;
    stream.opaque = Z_NULL;
//	stream.data_type = Z_TEXT;

    stream.next_in = compressBytes;
    stream.avail_in = cmpSize;
    /* Check for source > 64K on 16-bit machine: */
    if ((unsigned long)stream.avail_in != cmpSize) 
    {
		printf("Error: zlib_uncompress2: stream.avail_in != cmpSize");
		exit(1);
		return -1;
	}

    stream.next_out = *oriData;
    stream.avail_out = targetOriSize;
    //if ((uLong)stream.avail_out != *destLen) return Z_BUF_ERROR;

    int err = inflateInit(&stream);
    //int windowBits = 15;
    //int err = inflateInit2(&stream, windowBits);
    if (err != Z_OK)
    {
		printf("Error: zlib_uncompress2: err != Z_OK\n");
		return -1;
	}

    err = inflate(&stream, Z_FINISH);
    if (err != Z_STREAM_END) {
        inflateEnd(&stream);
        if (err == Z_NEED_DICT || (err == Z_BUF_ERROR && stream.avail_in == 0))
            return Z_DATA_ERROR;
        return err;
    }
    outSize = stream.total_out;
    inflateEnd(&stream);
    return outSize;
}



