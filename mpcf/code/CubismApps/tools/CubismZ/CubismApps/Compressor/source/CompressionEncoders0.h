/*
 *  CompressionEncoders0.h
 *  
 *
 *  Created by Panos Hadjidoukas on 01/15/17.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <zlib.h>	// always needed

inline int deflate_inplace(z_stream *strm, unsigned char *buf, unsigned len, unsigned *max);
inline size_t zdecompress(unsigned char * inputbuf, size_t ninputbytes, unsigned char * outputbuf, const size_t maxsize);


inline size_t zdecompress0(unsigned char * inputbuf, size_t ninputbytes, unsigned char * outputbuf, const size_t maxsize)
{
	int decompressedbytes = 0;

	decompressedbytes = ninputbytes;
	memcpy(outputbuf, inputbuf, ninputbytes);

	return decompressedbytes;
}

inline int deflate_inplace0(z_stream *strm, unsigned char *buf, unsigned len, unsigned *max)
{
	int compressedbytes = len;
	strm->total_out = compressedbytes;
	*max = compressedbytes;
        return Z_OK;
}

