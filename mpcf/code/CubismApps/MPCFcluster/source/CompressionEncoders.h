/*
 *  CompressionEncoders.h
 *  
 *
 *  Created by Diego Rossinelli on 3/27/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#ifndef _COMPRESSIONENCODERS_H_
#define _COMPRESSIONENCODERS_H_ 1

#pragma once

#include <zlib.h>
#if defined(_USE_LZ4_)
#include <lz4.h>
#endif
#if defined(_USE_FPZIP_)
extern "C"
{
#include "myfpzip.h"
}
#endif

inline int deflate_inplace(z_stream *strm, unsigned char *buf, unsigned len, unsigned *max);
inline size_t zdecompress(unsigned char * inputbuf, size_t ninputbytes, unsigned char * outputbuf, const size_t maxsize);


inline size_t zdecompress(unsigned char * inputbuf, size_t ninputbytes, unsigned char * outputbuf, const size_t maxsize)
{
	int decompressedbytes = 0;
#if defined(_USE_ZLIB_)
	z_stream datastream = {0};
	datastream.total_in = datastream.avail_in = ninputbytes;
	datastream.total_out = datastream.avail_out = maxsize;
	datastream.next_in = inputbuf;
	datastream.next_out = outputbuf;
	
	const int retval = inflateInit(&datastream);
	
	if (retval == Z_OK && inflate(&datastream, Z_FINISH))
	{
		decompressedbytes = datastream.total_out;
	}
	else
	{
		printf("ZLIB DECOMPRESSION FAILURE!!\n");
		abort();
	}
	
	inflateEnd(&datastream);
#elif defined(_USE_LZ4_)
	decompressedbytes = LZ4_uncompress_unknownOutputSize((char *)inputbuf, (char*) outputbuf, ninputbytes, maxsize);
	if (decompressedbytes < 0)
	{
		printf("LZ4 DECOMPRESSION FAILURE!!\n");
		abort();
	}
#else
	decompressedbytes = ninputbytes;
	memcpy(outputbuf, inputbuf, ninputbytes);
#endif
	return decompressedbytes;
}

/* THIS CODE SERVES US TO COMPRESS IN-PLACE. TAKEN FROM THE WEB
 * http://stackoverflow.com/questions/12398377/is-it-possible-to-have-zlib-read-from-and-write-to-the-same-memory-buffer
 * 
 * Compress buf[0..len-1] in place into buf[0..*max-1].  *max must be greater
 than or equal to len.  Return Z_OK on success, Z_BUF_ERROR if *max is not
 enough output space, Z_MEM_ERROR if there is not enough memory, or
 Z_STREAM_ERROR if *strm is corrupted (e.g. if it wasn't initialized or if it
 was inadvertently written over).  If Z_OK is returned, *max is set to the
 actual size of the output.  If Z_BUF_ERROR is returned, then *max is
 unchanged and buf[] is filled with *max bytes of uncompressed data (which is
 not all of it, but as much as would fit).
 
 Incompressible data will require more output space than len, so max should
 be sufficiently greater than len to handle that case in order to avoid a
 Z_BUF_ERROR. To assure that there is enough output space, max should be
 greater than or equal to the result of deflateBound(strm, len).
 
 strm is a deflate stream structure that has already been successfully
 initialized by deflateInit() or deflateInit2().  That structure can be
 reused across multiple calls to deflate_inplace().  This avoids unnecessary
 memory allocations and deallocations from the repeated use of deflateInit()
 and deflateEnd(). */
inline int deflate_inplace(z_stream *strm, unsigned char *buf, unsigned len,
						   unsigned *max)
{
#if defined(_USE_ZLIB_)
    int ret;                    /* return code from deflate functions */
    unsigned have;              /* number of bytes in temp[] */
    unsigned char *hold;        /* allocated buffer to hold input data */
    unsigned char temp[11];     /* must be large enough to hold zlib or gzip
								 header (if any) and one more byte -- 11
								 works for the worst case here, but if gzip
								 encoding is used and a deflateSetHeader()
								 call is inserted in this code after the
								 deflateReset(), then the 11 needs to be
								 increased to accomodate the resulting gzip
								 header size plus one */
	
	
    /* initialize deflate stream and point to the input data */
    ret = deflateReset(strm);
    
    if (ret != Z_OK)
        return ret;
    strm->next_in = buf;
    strm->avail_in = len;
	
    /* kick start the process with a temporary output buffer -- this allows
	 deflate to consume a large chunk of input data in order to make room for
	 output data there */
    if (*max < len)
        *max = len;
    strm->next_out = temp;
    strm->avail_out = sizeof(temp) > *max ? *max : sizeof(temp);
    ret = deflate(strm, Z_FINISH);
	
    if (ret == Z_STREAM_ERROR)
        return ret;
	
    /* if we can, copy the temporary output data to the consumed portion of the
	 input buffer, and then continue to write up to the start of the consumed
	 input for as long as possible */
    have = strm->next_out - temp;
    if (have <= (strm->avail_in ? len - strm->avail_in : *max)) {
        memcpy(buf, temp, have);
        strm->next_out = buf + have;
        have = 0;
        while (ret == Z_OK) {
            strm->avail_out = strm->avail_in ? strm->next_in - strm->next_out :
			(buf + *max) - strm->next_out;
            ret = deflate(strm, Z_FINISH);
        }
        if (ret != Z_BUF_ERROR || strm->avail_in == 0) {
            *max = strm->next_out - buf;
            return ret == Z_STREAM_END ? Z_OK : ret;
        }
    }
	
    /* the output caught up with the input due to insufficiently compressible
	 data -- copy the remaining input data into an allocated buffer and
	 complete the compression from there to the now empty input buffer (this
	 will only occur for long incompressible streams, more than ~20 MB for
	 the default deflate memLevel of 8, or when *max is too small and less
	 than the length of the header plus one byte) */
	
    hold = (unsigned char*)strm->zalloc(strm->opaque, strm->avail_in, 1);
    if (hold == Z_NULL)
        return Z_MEM_ERROR;
    memcpy(hold, strm->next_in, strm->avail_in);
    strm->next_in = hold;
    if (have) {
        memcpy(buf, temp, have);
        strm->next_out = buf + have;
    }
    strm->avail_out = (buf + *max) - strm->next_out;
    ret = deflate(strm, Z_FINISH);
    strm->zfree(strm->opaque, hold);
    *max = strm->next_out - buf;
    return ret == Z_OK ? Z_BUF_ERROR : (ret == Z_STREAM_END ? Z_OK : ret);

#elif defined(_USE_LZ4_)
	#define ZBUFSIZE (4*1024*1024)	/* fix this */
	static char bufzlib[ZBUFSIZE];	/* and this per thread (threadprivate or better an small cyclic array of buffers ) */

	if (ZBUFSIZE < *max) {
		printf("small ZBUFSIZE\n");
		abort();
	}

	int ninputbytes = len;          
	int compressedbytes;

#pragma omp critical
	{
	compressedbytes = LZ4_compress((char*) buf, (char *)bufzlib, ninputbytes);
	memcpy(buf, bufzlib, compressedbytes);	
	}

	strm->total_out = compressedbytes;
	*max = compressedbytes;
        return (compressedbytes >= 0) ? Z_OK : Z_BUF_ERROR;

#else
	int compressedbytes = len;
	strm->total_out = compressedbytes;
	*max = compressedbytes;
	return Z_OK;
#endif
}

#endif
