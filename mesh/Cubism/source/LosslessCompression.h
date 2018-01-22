/*
 *  LosslessCompression.h
 *  Cubism
 *
 *  Created by Panos Hadjidoukas on 3/20/14.
 *  Copyright 2014 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <cstring>
#include <zlib.h>	// always needed

#if defined(_USE_FPZIP_)
extern "C"
{
#include "myfpzip.h" 
}
#endif

inline size_t ZZcompress(unsigned char *buf, unsigned len, int layout[4], unsigned *max);
inline size_t ZZdecompress(unsigned char * inputbuf, size_t ninputbytes, int layout[4], unsigned char * outputbuf, const size_t maxsize);


inline size_t ZZdecompress(unsigned char * inputbuf, size_t ninputbytes, int layout[4], unsigned char * outputbuf, const size_t maxsize)
{
	int decompressedbytes = 0;

//	fpz_decompress4D((char *)inputbuf, ninputbytes, layout, (char *) outputbuf, (unsigned int *)&decompressedbytes, (sizeof(Real)==4)?1:0);
	fpz_decompress3D((char *)inputbuf, ninputbytes, layout, (char *) outputbuf, (unsigned int *)&decompressedbytes, (sizeof(Real)==4)?1:0);

#if DBG
	printf("fpz: %d to %d\n", ninputbytes, decompressedbytes);
#endif
	if (decompressedbytes < 0)
	{
		printf("FPZIP DECOMPRESSION FAILURE!!\n");
		abort();
	}
	return decompressedbytes;
}

inline size_t ZZcompress(unsigned char *buf, unsigned len, int layout[4], unsigned *max)
{
#if 1
	int zbufsize = 0;
	char *zbuf = NULL;
	zbufsize = len + 4096; 
	zbuf = (char *)malloc(zbufsize);
#else //keep the buffer allocated
	static int zbufsize = 0;
	static char *zbuf = NULL;
	if (zbuf == 0) {
		zbufsize = len + 4096; 
		zbuf = (char *)malloc(zbufsize);
	}
#endif

	if (zbufsize < *max) {
		printf("small ZBUFSIZE\n");
		abort();
	}

	int ninputbytes = len;          
	int compressedbytes;

//	fpz_compress4D((void *) buf, ninputbytes, layout, (void *) zbuf, (unsigned int *)&compressedbytes, (sizeof(Real)==4)?1:0);
	fpz_compress3D((void *) buf, ninputbytes, layout, (void *) zbuf, (unsigned int *)&compressedbytes, (sizeof(Real)==4)?1:0);
	memcpy(buf, zbuf, compressedbytes);	

	*max = compressedbytes;
#if 1
	free(zbuf);
#endif
	return compressedbytes;
}
