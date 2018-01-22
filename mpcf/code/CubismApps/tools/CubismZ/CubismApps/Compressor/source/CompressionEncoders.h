/*
 *  CompressionEncoders.h
 *
 *
 *  Created by Diego Rossinelli on 3/27/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#if defined(_USE_ZSTD_)
#include "zstd_zlibwrapper.h"
#else
#include <zlib.h>	// always needed
#endif

#if defined(_USE_LZ4_)
#include <lz4.h>
#endif

#if defined(_USE_LZF_)
extern "C"
{
#include "lzf.h"
}
#endif

#if defined(_USE_LZMA_)
extern "C"
{
#include "mylzma.h"
}
#endif

#if defined(_USE_ZSTD0_)
extern "C"
{
#include "zstd.h"
}
#endif

#if defined(_USE_BLOSC_)
#include "blosc.h"
#endif

#if defined(_USE_ZOPFLI_)
extern "C"
{
#include "myzopfli.h"
}
#endif

#if defined(_USE_FPC_)||defined(_USE_FPC2_)
extern "C"
{
#include "fpc.h"
}
#endif

#if defined(_USE_FPZIP_)||defined(_USE_FPZIP2_)
extern "C"
{
#include "myfpzip.h"
}
#endif

#if defined(_USE_DRAIN_)
#include "drain.h"
#endif

#if defined(_USE_ZFP_)||defined(_USE_ZFP3_)
#include "myzfp.h"
#endif

#if defined(_USE_SZ_)||defined(_USE_SZ3_)
extern "C"
{
#include "rw.h"
#include "sz.h"
}
#endif

#if defined(_USE_ISA_)
#include "myisa.h"
#endif

inline int deflate_inplace(z_stream *strm, unsigned char *buf, unsigned len, unsigned *max);
inline size_t zdecompress(unsigned char * inputbuf, size_t ninputbytes, unsigned char * outputbuf, const size_t maxsize);


inline size_t zdecompress(unsigned char * inputbuf, size_t ninputbytes, unsigned char * outputbuf, const size_t maxsize)
{
#if defined(VERBOSE)
	printf("zdecompress has been called for %d input bytes\n", ninputbytes);
#endif

	int decompressedbytes = 0;
#if defined(_USE_ZLIB_) || defined(_USE_ZSTD_)
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
#elif defined(_USE_LZF_)
	decompressedbytes = lzf_decompress((const void *)inputbuf, ninputbytes, (void *) outputbuf, maxsize);
	if (decompressedbytes < 0)
	{
		printf("LZF DECOMPRESSION FAILURE!!\n");
		abort();
	}
#elif defined(_USE_LZMA_)
#if 0
	unsigned char *lzma_decompressed;
	int rc = simpleDecompress(ELZMA_lzma, inputbuf, ninputbytes, &lzma_decompressed, (size_t *)&decompressedbytes);
	if (rc != ELZMA_E_OK)
	{
		printf("LZMA DECOMPRESSION FAILURE!!\n");
		abort();
	}
	memcpy(outputbuf, lzma_decompressed, decompressedbytes);
	free(lzma_decompressed);
#else
	decompressedbytes = lzma_decompress((unsigned char *) inputbuf, ninputbytes, (unsigned char *) outputbuf, maxsize);
	if (decompressedbytes < 0)
	{
		printf("LZMA DECOMPRESSION FAILURE!!\n");
		abort();
	}
#endif
#elif defined(_USE_ZSTD0_)
	{
	ZSTD_DStream* const dstream = ZSTD_createDStream();
	if (dstream==NULL) { fprintf(stderr, "ZSTD_createDStream() error \n"); exit(10); }
	size_t const initResult = ZSTD_initDStream(dstream);
	if (ZSTD_isError(initResult)) { fprintf(stderr, "ZSTD_initDStream() error : %s \n", ZSTD_getErrorName(initResult)); exit(11); }

	ZSTD_inBuffer input = { inputbuf, ninputbytes, 0};
	ZSTD_outBuffer output = { outputbuf, maxsize, 0};
	size_t toRead = ZSTD_decompressStream(dstream, &output, &input);  /* toRead : size of next compressed block */
	if (ZSTD_isError(toRead)) { fprintf(stderr, "ZSTD_decompressStream() error : %s \n", ZSTD_getErrorName(toRead)); exit(12); }
	decompressedbytes = output.pos;

	ZSTD_freeDStream(dstream);
	//printf("zstd0: %ld -> %ld  (%.2lfx)\n", ninputbytes, decompressedbytes, (1.0*decompressedbytes)/ninputbytes);
	}
#elif defined(_USE_BLOSC_)
	/* Decompress  */
	size_t dsize = blosc_decompress(inputbuf, outputbuf, maxsize);
	if (dsize < 0) {
		printf("Decompression error.  Error code: %d\n", dsize);
	}
	decompressedbytes = dsize;
#elif defined(_USE_ZOPFLI_)
	decompressedbytes = zopfli_decompress((const void *) inputbuf, ninputbytes, (void *) outputbuf, maxsize);
	if (decompressedbytes < 0)
	{
		printf("ZOPFLI DECOMPRESSION FAILURE!!\n");
		abort();
	}
#elif defined(_USE_FPC2_)
	fpc_decompress((char *)inputbuf, ninputbytes, (char*) outputbuf, &decompressedbytes);
	if (decompressedbytes < 0)
	{
		printf("FPC DECOMPRESSION FAILURE!!\n");
		abort();
	}
#elif defined(_USE_FPZIP2_)	/* does not work as the buffer contains a mix of integers + floats */
	fpz_decompress1D((char *)inputbuf, ninputbytes, (char *) outputbuf, (unsigned int *)&decompressedbytes, (sizeof(Real)==4)?1:0);
	printf("fpz: %d to %d\n", ninputbytes, decompressedbytes);
	if (decompressedbytes < 0)
	{
		printf("FPZIP DECOMPRESSION FAILURE!!\n");
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
#if defined(_USE_ZLIB_)||defined(_USE_ZSTD_)
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

#elif defined(_USE_LZ4_)||defined(_USE_LZF_)||defined(_USE_LZMA_)||defined(_USE_ZOPFLI_)||defined(_USE_FPC2_)||defined(_USE_FPZIP2_)||defined(_USE_ZSTD0_)||defined(_USE_BLOSC_)

	#define ZBUFSIZE (4*1024*1024)	/* fix this */

//	static char bufzlib[ZBUFSIZE];	/* and this per thread (threadprivate or better a small cyclic array of buffers ) */
	#define MAXBUFFERS	12	/* todo */
	static char bufzlibA[MAXBUFFERS][ZBUFSIZE];	/* and this per thread (threadprivate or better a small cyclic array of buffers ) */

	if (omp_get_num_threads() > MAXBUFFERS) {
		printf("small number of buffers (%d)\n", MAXBUFFERS);
		abort();
	}

	if (ZBUFSIZE < *max) {
		printf("small ZBUFSIZE\n");
		abort();
	}

	char *bufzlib = (char *)&bufzlibA[omp_get_thread_num()][0];

	int ninputbytes = len;
	int compressedbytes;

#if (MAXBUFFERS == 1)
#pragma omp critical
#endif
	{
#if defined(_USE_LZ4_)
		compressedbytes = LZ4_compress((char*) buf, (char *) bufzlib, ninputbytes);
#elif defined(_USE_LZF_)
		compressedbytes = lzf_compress((void *) buf, ninputbytes, (void *) bufzlib, ZBUFSIZE);
#elif defined(_USE_LZMA_)
#if 0
		int rc;
		unsigned char *lzma_compressed;
		rc = simpleCompress(ELZMA_lzma, (unsigned char *) buf, ninputbytes, &lzma_compressed, (size_t *)&compressedbytes);
//		printf("lzma compression: %ld -> %ld\n", ninputbytes, compressedbytes);
		if (rc != ELZMA_E_OK) {
			compressedbytes = -1;
//			printf("rz != ELZMA_E_OK\n"); abort();
		}
		memcpy(bufzlib, lzma_compressed, compressedbytes);
		free(lzma_compressed);
#else
		compressedbytes = lzma_compress((unsigned char *) buf, ninputbytes, (unsigned char *)bufzlib, ZBUFSIZE);
#endif
#elif defined(_USE_ZSTD0_)
		{
		int cLevel = 1;
		ZSTD_CStream* const cstream = ZSTD_createCStream();
		size_t const initResult = ZSTD_initCStream(cstream, cLevel);

		ZSTD_inBuffer input = { buf, ninputbytes, 0 };
		ZSTD_outBuffer output = { bufzlib, ZBUFSIZE, 0 };
		size_t toRead = ZSTD_compressStream(cstream, &output, &input);
		if (ZSTD_isError(toRead)) { fprintf(stderr, "ZSTD_compressStream() error : %s \n", ZSTD_getErrorName(toRead)); exit(12); }

		//ZSTD_outBuffer output = { buffOut, buffOutSize, 0 };
		size_t const remainingToFlush = ZSTD_endStream(cstream, &output);   /* close frame */
		if (remainingToFlush) { fprintf(stderr, "not fully flushed"); exit(13); }

		compressedbytes = output.pos;
		ZSTD_freeCStream(cstream);

		//printf("zstd0: %ld -> %ld  (%.2lfx)\n", ninputbytes, compressedbytes, (1.0*ninputbytes)/compressedbytes);

		}
#elif defined(_USE_BLOSC_)
		{
		int clevel = 5;
		int doshuffle = 1;
		int typesize = sizeof(Real);

		/* Compress with clevel=5 and shuffle active  */
		size_t csize = blosc_compress(clevel, doshuffle, typesize, ninputbytes, buf, bufzlib, ZBUFSIZE);
		if (csize < 0) {
			printf("Compression error.  Error code: %d\n", csize);
		}
		compressedbytes = csize;
		//printf("blosc: %ld -> %ld\n", ninputbytes, compressbytes);
		}
#elif defined(_USE_ZOPFLI_)
		compressedbytes = zopfli_compress((const unsigned char *) buf, ninputbytes, (unsigned char *) bufzlib);
#elif defined(_USE_FPC2_)
		fpc_compress((char *) buf, ninputbytes, (char *) bufzlib, &compressedbytes, 10);
#elif defined(_USE_FPZIP2_)
		fpz_compress1D((void *) buf, ninputbytes, (void *) bufzlib, (unsigned int *)&compressedbytes, (sizeof(Real)==4)?1:0);
#endif
		memcpy(buf, bufzlib, compressedbytes);
	}

	strm->total_out = compressedbytes;
	*max = compressedbytes;
        return (compressedbytes >= 0) ? Z_OK : Z_BUF_ERROR;

#else	/* NO COMPRESSION */

	int compressedbytes = len;
	strm->total_out = compressedbytes;
	*max = compressedbytes;
        return Z_OK;
#endif
}


#if defined(_USE_SHUFFLE_)
void shuffle(char *in, int n, int s)
{
        int i, j, k;
        int b;

        char *tmp = (char *)malloc(n);

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

        char *tmp = (char *)malloc(n);
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
#endif
