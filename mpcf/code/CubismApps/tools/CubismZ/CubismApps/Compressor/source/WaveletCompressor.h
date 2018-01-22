/*
 *  WaveletCompressor.h
 *
 *
 *  Created by Diego Rossinelli on 3/4/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#ifndef _WAVELETCOMPRESSOR_H_
#define _WAVELETCOMPRESSOR_H_ 1

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cstdio>

#include "FullWaveletTransform.h"

#if defined(_USE_ZSTD_)
#include "zstd_zlibwrapper.h"
#else
#include <zlib.h>	// always needed
#endif

#if defined(_USE_LZ4_)
#include <lz4.h>
#endif

#ifdef _SP_COMP_
typedef float Real;
#else
typedef double Real;
#endif

static const bool lifting_scheme = false; // peh: disabled: xxx: check this again

template<int DATASIZE1D, typename DataType>
class WaveletCompressorGeneric
{
protected:

	enum
	{
		BS3 = DATASIZE1D * DATASIZE1D * DATASIZE1D,
		BITSETSIZE = (BS3 + 7) / 8,
//		BITSETSIZE = BS3,
		BUFMAXSIZE = BITSETSIZE + sizeof(DataType) * BS3
	};

	WaveletsOnInterval::FullTransform<DATASIZE1D> full;

private:

	unsigned char bufcompression[BUFMAXSIZE];

	size_t bufsize;

public:

	WaveletsOnInterval::FwtAp (& uncompressed_data()) [DATASIZE1D][DATASIZE1D][DATASIZE1D] { return full.data; }

	virtual void * compressed_data() { return bufcompression; }

	virtual	size_t compress(const float threshold, const bool float16, int wtype);
	virtual	size_t compress(const float threshold, const bool float16, bool swap, int wtype);

	virtual void decompress(const bool float16, size_t bytes, int wtype);

	virtual void decompress(const bool float16, size_t ninputbytes, int wtype, DataType data[DATASIZE1D][DATASIZE1D][DATASIZE1D])
	{
		decompress(float16, ninputbytes, wtype);

		this->copy_to(data);
	}

	void copy_to(DataType data[DATASIZE1D][DATASIZE1D][DATASIZE1D])
	{
		DataType * const dst = &data[0][0][0];
		const WaveletsOnInterval::FwtAp * const src = &full.data[0][0][0];
		for(int i = 0; i < BS3; ++i)
		{
			dst[i] = src[i];
			assert(!std::isnan(dst[i]));
		}
	}

	void copy_from(const DataType data[DATASIZE1D][DATASIZE1D][DATASIZE1D])
	{
		const DataType * const src = &data[0][0][0];
		WaveletsOnInterval::FwtAp * const dst = &full.data[0][0][0];
		for(int i = 0; i < BS3; ++i)
		{
			dst[i] = src[i];
			assert(!std::isnan(dst[i]));
		}
	}
};

template<int DATASIZE1D, typename DataType>
class WaveletCompressorGeneric_zlib: public WaveletCompressorGeneric<DATASIZE1D, DataType>
{
	unsigned char bufzlib[WaveletCompressorGeneric<DATASIZE1D, DataType>::BUFMAXSIZE];

public:

	void * compressed_data() { return bufzlib; }

	size_t compress(const float threshold, const bool float16, int wtype)
	{
		int compressedbytes = 0;

		const size_t ninputbytes = WaveletCompressorGeneric<DATASIZE1D, DataType>::compress(threshold, float16, wtype);

#if defined(_USE_ZLIB_) || defined(_USE_ZSTD_)
		z_stream datastream = { 0 };
		datastream.total_in = datastream.total_out = 0;
		datastream.avail_in = ninputbytes;
		datastream.avail_out = WaveletCompressorGeneric<DATASIZE1D, DataType>::BUFMAXSIZE;
		datastream.next_in = (unsigned char*) WaveletCompressorGeneric<DATASIZE1D, DataType>::compressed_data();
		datastream.next_out = bufzlib;

//		if (Z_OK == deflateInit(&datastream, Z_DEFAULT_COMPRESSION) && Z_STREAM_END == deflate(&datastream, Z_FINISH))
		if (Z_OK == deflateInit(&datastream, Z_BEST_COMPRESSION) && Z_STREAM_END == deflate(&datastream, Z_FINISH))
			compressedbytes = datastream.total_out;
		else
		{
			printf("ZLIB COMPRESSION FAILURE!!\n");
			abort();
		}

		deflateEnd( &datastream );
#elif defined(_USE_LZ4)
		compressedbytes = LZ4_compress((char*) WaveletCompressorGeneric<DATASIZE1D, DataType>::compressed_data(), (char *)bufzlib, ninputbytes);
		if (compressedbytes < 0)
		{
			printf("LZ4 COMPRESSION FAILURE!!\n");
			abort();
		}
#endif

		return compressedbytes;
	}

	size_t compress(const float threshold, const bool float16, bool swap, int wtype)
	{
		int compressedbytes = 0;

		const size_t ninputbytes = WaveletCompressorGeneric<DATASIZE1D, DataType>::compress(threshold, float16, swap, wtype);

#if defined(_USE_ZLIB_) || defined(_USE_ZSTD_)
		z_stream datastream = { 0 };
		datastream.total_in = datastream.total_out = 0;
		datastream.avail_in = ninputbytes;
		datastream.avail_out = WaveletCompressorGeneric<DATASIZE1D, DataType>::BUFMAXSIZE;
		datastream.next_in = (unsigned char*) WaveletCompressorGeneric<DATASIZE1D, DataType>::compressed_data();
		datastream.next_out = bufzlib;

//		if (Z_OK == deflateInit(&datastream, Z_DEFAULT_COMPRESSION) && Z_STREAM_END == deflate(&datastream, Z_FINISH))
		if (Z_OK == deflateInit(&datastream, Z_BEST_COMPRESSION) && Z_STREAM_END == deflate(&datastream, Z_FINISH))
			compressedbytes = datastream.total_out;
		else
		{
			printf("ZLIB COMPRESSION FAILURE!!\n");
			abort();
		}

		deflateEnd( &datastream );
#elif defined(_USE_LZ4)
		compressedbytes = LZ4_compress((char*) WaveletCompressorGeneric<DATASIZE1D, DataType>::compressed_data(), (char *)bufzlib, ninputbytes);
		if (compressedbytes < 0)
		{
			printf("LZ4 COMPRESSION FAILURE!!\n");
			abort();
		}
#endif

		return compressedbytes;
	}

	void decompress(const bool float16, size_t ninputbytes, int wtype)
	{
		int decompressedbytes = 0;

#if defined(_USE_ZLIB_) || defined(_USE_ZSTD_)
		z_stream datastream = { 0 };
		datastream.total_in = datastream.total_out = 0;
		datastream.avail_in = ninputbytes;
		datastream.avail_out = WaveletCompressorGeneric<DATASIZE1D, DataType>::BUFMAXSIZE;
		datastream.next_in = bufzlib;
		datastream.next_out = (unsigned char*) WaveletCompressorGeneric<DATASIZE1D, DataType>::compressed_data();

		const int retval = inflateInit(&datastream);

		if (retval == Z_OK && inflate(&datastream, Z_FINISH))
			decompressedbytes = datastream.total_out;
                else
                {
                        printf("ZLIB DECOMPRESSION FAILURE!!\n");
                        abort();
                }

		inflateEnd(&datastream);

		WaveletCompressorGeneric<DATASIZE1D, DataType>::decompress(float16, decompressedbytes, wtype);
#elif defined(_USE_LZ4)
                decompressedbytes = LZ4_uncompress_unknownOutputSize((char *)bufzlib, (char*) WaveletCompressorGeneric<DATASIZE1D, DataType>::compressed_data(),
									ninputbytes, WaveletCompressorGeneric<DATASIZE1D, DataType>::BUFMAXSIZE);
		if (decompressedbytes < 0)
		{
			printf("LZ4 DECOMPRESSION FAILURE!!\n");
			abort();
		}
#endif
	}
};

#ifdef _BLOCKSIZE_
typedef WaveletCompressorGeneric<_BLOCKSIZE_, Real> WaveletCompressor;
typedef WaveletCompressorGeneric_zlib<_BLOCKSIZE_, Real> WaveletCompressor_zlib;
#endif


#endif
