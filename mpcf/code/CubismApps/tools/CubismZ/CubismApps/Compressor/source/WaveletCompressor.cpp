/*
 *  WaveletCompressor.cpp
 *  
 *
 *  Created by Diego Rossinelli on 3/4/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#ifdef _SP_COMP_
typedef float Real;
#else
typedef double Real;
#endif

#include <cstdio>
#include <bitset>
#include <cassert>

using namespace std;

#include "WaveletCompressor.h"

//#define _USE_FLOAT16_	// peh: put it as an option in Makefile.config


#if defined(_USE_SPDP3_)
#include "myspdp.h"
//static size_t spdp_compress_data( const char* buffIn, size_t buffInSize, char* buffOut, size_t buffOutSize )
//static size_t spdp_uncompress_data( const char* buffIn, size_t buffInSize, char* buffOut, size_t buffOutSize)
#endif

#if defined(_USE_FPZIP3_)
#include "myfpzip.h"
#endif

#if defined(_USE_ZFP3_)
#include "myzfp.h"
#endif

#if defined(_USE_SZ3_)
extern "C"
{
#include "rw.h"
#include "sz.h"
}
#endif


void swapbytes(unsigned char *mem, int nbytes)
{
        unsigned char buf[8];
        for (int i = 0; i < nbytes; i++) buf[i] = mem[i];
        for (int i = 0; i < nbytes; i++) mem[nbytes-i-1] = buf[i];
}



#if defined(_USE_ZEROBITS_)
void float_zero_bits3(unsigned int *ul, int n)
{
        unsigned int _ul = *ul;
        unsigned long mask = ~((1 << n) - 1);
        _ul = _ul & mask;

#if 0
        if (n == 0)
                _ul = _ul;
        else if (n == 4)
                _ul = _ul & 0xfffffff0;
        else if (n == 8)
                _ul = _ul & 0xffffff00;
        else if (n == 12)
                _ul = _ul & 0xfffff000;
        else if (n == 16)
                _ul = _ul & 0xffff0000;
#endif

        *ul = _ul;


}

float float_zero_bits3b(float x, int bits)
{
        union { int ix; float fx; } v;
        v.fx = x;
        v.ix &= ~((1 << bits) - 1);

        return v.fx;
}
#endif

#if defined(_USE_SHUFFLE3_)
static void shuffle3(char *in, int n, int s)
{
	int i, j, b;
	char *tmp = (char *)malloc(n);

	j = 0;
	for (b = 0; b < s; b++)
	{
		for (i = b; i < n; i+= s)
		{
			tmp[j] = in[i];
			j++;
		}
	}

	memcpy(in, tmp, n);
	free(tmp);
}

void reshuffle3(char *in, int n, int s)
{
	int i, j, b;
	char *tmp = (char *)malloc(n);

	int c = n/s;

	j = 0;
	for (i = 0; i < n/s; i++)
	{
		for (b = 0; b < s; b++)
		{
			tmp[j] = in[i + b*c]; j++;
		}
	}

	memcpy(in, tmp, n);
	free(tmp);
}

#endif

template<int N>
void serialize_bitset(bitset<N> mybits, unsigned char * const buf, const int nbytes)
{
	assert(nbytes == (N + 7) / 8);
	
	const int nicebits = 8 * (N / 8);
	
	for(int i = 0, B = 0; i < nicebits; i += 8, ++B)
	{
		unsigned char c = 0;
		
		for(int b = 0; b < 8; ++b)
			c |= mybits[i + b] << b;
		
		buf[B] = c;
	}
	
	if (nicebits < N)
	{
		unsigned char c = 0;
		
		for(int b = nicebits; b < N; ++b)
			c |= mybits[b] << (b - nicebits);
		
		buf[nbytes - 1] = c;
	}
}

template<int N>
int deserialize_bitset(bitset<N>& mybits, const unsigned char * const buf, const int nbytes)
{
	assert(nbytes == (N + 7) / 8);
	
	int sum = 0;
	
	const int nicebits = 8 * (N / 8);
	
	for(int i = 0, B = 0; i < nicebits; i += 8, ++B)
	{
		const unsigned char c = buf[B];
		
		for(int b = 0; b < 8; ++b)
			sum += (mybits[i + b] = (c >> b) & 1);		
	}
	
	if (nicebits < N)
	{
		const unsigned char c = buf[nbytes - 1];
		
		for(int b = nicebits; b < N; ++b)
			sum += (mybits[b] = (c >> (b - nicebits)) & 1);
	}
	
	return sum;
}

unsigned short _cvt2f16(const float xfloat)
{
	assert(sizeof(unsigned short) == 2);
	
	const unsigned int x = *(unsigned int *)& xfloat;
	
	//copy sign
	unsigned short retval = (x & 0x80000000) >> 16; 
	
	//reshape exponent
	const int e = ((x & 0x7fffffff) >> 23) - 127;
	retval |= max(0, min(31, e + 15)) << 10;
	
	//reshape mantissa
	const unsigned int m = x & 0x007fffff;
	retval |= m >> 13;
	
	return retval;
}

float _cvtfromf16(const unsigned short f16)
{
	//copy sign
	int retval = (f16 & 0x8000) << 16;
	
	//reshape exponent
	const int e = (((0x1f << 10) & f16) >> 10) - 15;
	retval |= (e + 127) << 23;
	
	//reshape mantissa
	const int m = f16 & 0x3ff;
	retval |= m << 13;
	
	return *(float *)& retval;
}

template<int DATASIZE1D, typename DataType>
size_t WaveletCompressorGeneric<DATASIZE1D, DataType>::compress(const float threshold, const bool float16, int wtype)//, const DataType data[DATASIZE1D][DATASIZE1D][DATASIZE1D])
{				
	full.fwt(wtype);
	
	assert(BITSETSIZE % sizeof(DataType) == 0);
	
	bitset<BS3> mask;
	const int survivors = full.template threshold<DataType, DATASIZE1D>(threshold, mask, (DataType *)(bufcompression + BITSETSIZE));

#if 0
	//compare it against the old code
	{
		pair<vector<DataType> , bitset<BS3> > _survivors = full.template threshold<DataType>(threshold);
		assert(mask == _survivors.second);
		vector<DataType> newvals( (DataType *)(bufcompression + BITSETSIZE),  survivors + (DataType *)(bufcompression + BITSETSIZE));
		
		assert(newvals == _survivors.first);
	}
#endif

	serialize_bitset<BS3>(mask, bufcompression, BITSETSIZE);

#if defined(_USE_SPDP3_)||defined(_USE_FPZIP3_)||defined(_USE_ZFP3_)||defined(_USE_SZ3_)
	{
	int survbytes = sizeof(DataType)*survivors;
	int outbytes;
	int bufsize = 4*1024*1024*sizeof(char);
	unsigned char *tmp = (unsigned char *)malloc(bufsize);

#if defined(_USE_SPDP3_)

#if 0
	float *s = (float *) (bufcompression + BITSETSIZE);
	printf("s = [\n");
	for (int i = 0; i < survivors; i++) printf("%f\n", s[i]);
	printf("]\n");
#endif

	outbytes = spdp_compress_data((char *)(bufcompression + BITSETSIZE), survbytes, (char *)tmp, bufsize);
#if 0
	printf("spdp3 (%d): %d -> %d (%.2fx)\n", survivors, survbytes, outbytes, 1.0*survbytes/outbytes);
#endif

#if 0
	printf("buf = {\n");
	for (int i = 0; i < outbytes; i++) printf("%d,\n", tmp[i]);
	printf("}\n");
#endif

#elif defined(_USE_FPZIP3_)
	int layout[4] = {1, 1, survivors, 1};
	int is_float = sizeof(Real)==4;

	int fpzip_bits = 8*sizeof(Real);
	if(getenv("FPZIP_BITS")) fpzip_bits = atoi(getenv("FPZIP_BITS"));
	if (fpzip_bits > 8*sizeof(Real)) fpzip_bits = 8*sizeof(Real);

//	fpz_compress3D((char *) (bufcompression + BITSETSIZE), survbytes, layout, tmp, (unsigned int *)&outbytes, 1, 8*sizeof(Real));
	fpz_compress3D((char *) (bufcompression + BITSETSIZE), survbytes, layout, tmp, (unsigned int *)&outbytes, is_float, fpzip_bits);
#if VERBOSE
	printf("fpz3 : %d -> %d (%.2fx)\n", survbytes, outbytes, 1.0*survbytes/outbytes);
#endif

#elif defined(_USE_ZFP3_)
	double zfp_acc = 0.01;
	if(getenv("ZFP_ACC")) zfp_acc = atof(getenv("ZFP_ACC"));
	int is_float = sizeof(Real)==4;
	int layout[4] = {survivors, 1, 1, 1};
	size_t nbytes_zfp;
	int status = zfp_compress_buffer((char *)(bufcompression + BITSETSIZE), layout[0], layout[1], layout[2], zfp_acc, is_float, (unsigned char *)tmp, &nbytes_zfp);
	outbytes = nbytes_zfp;
#if VERBOSE
	printf("zfp3 : %d -> %d (%.2fx)\n", survbytes, outbytes, 1.0*survbytes/outbytes);
#endif
	if (outbytes > survbytes) {
		printf("zfp3 : %d -> %d (%.2fx)\n", survbytes, outbytes, 1.0*survbytes/outbytes);
		exit(1);
	}

#elif defined(_USE_SZ3_)
	double sz_abs_acc = 0.0;
	double sz_rel_acc = 0.0;
	if(getenv("SZ_ABS_ACC")) sz_abs_acc = atof(getenv("SZ_ABS_ACC"));
	//SZ_ABS_ACC=$PARAM
	//sz_abs_acc = (double) this->threshold;
	int layout[5] = {survivors, 1, 1, 0, 0};
	int bytes_sz;
	unsigned char *compressed_sz = SZ_compress_args(SZ_FLOAT, (unsigned char *)((char *)bufcompression + BITSETSIZE), &bytes_sz, ABS, sz_abs_acc, sz_rel_acc, layout[4], layout[3], layout[2], layout[1], layout[0]);
	outbytes = bytes_sz;
	memcpy(tmp, compressed_sz, outbytes);
	free(compressed_sz);
	printf("sz3 : %d -> %d (%.2fx)\n", survbytes, outbytes, 1.0*survbytes/outbytes);
#endif

	memcpy((char *) (bufcompression + BITSETSIZE), tmp, outbytes);
	free(tmp);

	return (BITSETSIZE + outbytes);
	}

#elif defined(_USE_SHUFFLE3_)||defined(_USE_ZEROBITS_)

 #if !defined(_USE_FLOAT16_)	// ZEROBITS and/or SHUFFLE

  #if defined(_USE_ZEROBITS_)
	// set some bits to zero
	for(int i = 0; i < survivors; ++i) 
	{
		DataType *ps = (i + (DataType *)(bufcompression + BITSETSIZE));
		float_zero_bits3((unsigned int *)ps, _ZEROBITS_);	// xxx: extend it for doubles
	}
  #endif

  #if defined(_USE_SHUFFLE3_)
	shuffle3((char *)(bufcompression + BITSETSIZE), survivors*sizeof(DataType), sizeof(DataType));
  #endif

 #else	// _USE_FLOAT16_ (+shuffle3)

	for(int i = 0; i < survivors; ++i) 
		//dangerous, but it looks like i know where i am going with this
		*(i + (unsigned short *)(bufcompression + BITSETSIZE)) = 
		_cvt2f16(*(i + (DataType *)(bufcompression + BITSETSIZE)));

	shuffle3((char *)(bufcompression + BITSETSIZE), survivors*sizeof(unsigned short), sizeof(unsigned short));

	return BITSETSIZE + sizeof(unsigned short) * survivors;

 #endif

#endif

	if (!float16) return BITSETSIZE + sizeof(DataType) * survivors;
	
	for(int i = 0; i < survivors; ++i) 
		//dangerous, but it looks like i know where i am going with this
		*(i + (unsigned short *)(bufcompression + BITSETSIZE)) = 
		_cvt2f16(*(i + (DataType *)(bufcompression + BITSETSIZE)));
	
	return BITSETSIZE + sizeof(unsigned short) * survivors;
}


template<int DATASIZE1D, typename DataType>
size_t WaveletCompressorGeneric<DATASIZE1D, DataType>::compress(const float threshold, const bool float16, bool swap, int wtype)//, const DataType data[DATASIZE1D][DATASIZE1D][DATASIZE1D])
{				
	full.fwt(wtype);
	
	exit(1);
	assert(BITSETSIZE % sizeof(DataType) == 0);
	
	bitset<BS3> mask;
	const int survivors = full.template threshold<DataType, DATASIZE1D>(threshold, mask, (DataType *)(bufcompression + BITSETSIZE));

#if 0
	//compare it against the old code
	{
		pair<vector<DataType> , bitset<BS3> > _survivors = full.template threshold<DataType>(threshold);
		assert(mask == _survivors.second);
		vector<DataType> newvals( (DataType *)(bufcompression + BITSETSIZE),  survivors + (DataType *)(bufcompression + BITSETSIZE));
		
		assert(newvals == _survivors.first);
	}
#endif

	serialize_bitset<BS3>(mask, bufcompression, BITSETSIZE);

#if defined(_USE_SHUFFLE3_)||defined(_USE_ZEROBITS_)

 #if !defined(_USE_FLOAT16_)	// ZEROBITS and/or SHUFFLE

  #if defined(_USE_ZEROBITS_)
	// set some bits to zero
	for(int i = 0; i < survivors; ++i) 
	{
		DataType *ps = (i + (DataType *)(bufcompression + BITSETSIZE));
		float_zero_bits3((unsigned int *)ps, _ZEROBITS_);	// xxx: extend it for doubles
	}
  #endif


	// peh: need to verify if this is the correct place for byte swapping. However, swapping is in general useless for the compression tool.
	if (swap)
	{
		unsigned char *buf = ((unsigned char *)bufcompression + BITSETSIZE);
		int sz = sizeof(DataType);	// = 4
		for(int i = 0; i < survivors; ++i) {
			swapbytes((unsigned char *)buf+i*sz, sz);
		}
	}


  #if defined(_USE_SHUFFLE3_)
	shuffle3((char *)(bufcompression + BITSETSIZE), survivors*sizeof(DataType), sizeof(DataType));
  #endif

 #else	// _USE_FLOAT16_ (+shuffle3)

	for(int i = 0; i < survivors; ++i) 
		//dangerous, but it looks like i know where i am going with this
		*(i + (unsigned short *)(bufcompression + BITSETSIZE)) = 
		_cvt2f16(*(i + (DataType *)(bufcompression + BITSETSIZE)));


	// peh: need to verify if this is the correct place for byte swapping. However, swapping is in general useless for the compression tool.
	if (swap)
	{
		unsigned char *buf = ((unsigned char *)bufcompression + BITSETSIZE);
		int sz = sizeof(DataType);	// = 4
		for(int i = 0; i < survivors; ++i) {
			swapbytes((unsigned char *)buf+i*sz, sz);
		}
	}

	shuffle3((char *)(bufcompression + BITSETSIZE), survivors*sizeof(unsigned short), sizeof(unsigned short));

	return BITSETSIZE + sizeof(unsigned short) * survivors;

 #endif

#endif

	if (!float16) return BITSETSIZE + sizeof(DataType) * survivors;
	
	for(int i = 0; i < survivors; ++i) 
		//dangerous, but it looks like i know where i am going with this
		*(i + (unsigned short *)(bufcompression + BITSETSIZE)) = 
		_cvt2f16(*(i + (DataType *)(bufcompression + BITSETSIZE)));
	
	return BITSETSIZE + sizeof(unsigned short) * survivors;
}


template<int DATASIZE1D, typename DataType>
void WaveletCompressorGeneric<DATASIZE1D, DataType>::decompress(const bool float16, size_t bytes, int wtype)//, DataType data[DATASIZE1D][DATASIZE1D][DATASIZE1D])
{
	bitset<BS3> mask;
	const int expected = deserialize_bitset<BS3>(mask, bufcompression, BITSETSIZE);

#if defined(_USE_SPDP3_)||defined(_USE_FPZIP3_)||defined(_USE_ZFP3_)||defined(_USE_SZ3_)
        {
	int outbytes;
	int bufsize = 4*1024*1024*sizeof(char);
	int survbytes = bytes - BITSETSIZE;
	int survivors = expected; //survbytes/sizeof(Real);

	//printf("compressed survbytes = %d, expected = %d\n", survbytes, expected);

	char *tmp = (char *)malloc(bufsize);
#if defined(_USE_SPDP3_)

#if 0
	unsigned char *tmp0 = bufcompression + BITSETSIZE;
	printf("buf = [\n");
	for (int i = 0; i < survbytes; i++) printf("%d\n", tmp0[i]);
	printf("]\n");
#endif

	outbytes = spdp_uncompress_data((char *)bufcompression + BITSETSIZE, bytes - BITSETSIZE, tmp, bufsize);
	//printf("spdp3(D) : %d -> %d\n", survbytes, outbytes);

	if (0)
	{
	float *s1 = (float *)tmp;
	for (int i = 0; i < survivors; i++)
	{
		printf("s1[%d] = %f\n", i, s1[i]);
	}
	}

#elif defined(_USE_FPZIP3_)
	//fpz_decompress3D(char *in, unsigned int inbytes, int layout[4], char *out, unsigned int *outbytes, int isfloat, int prec)
	int layout[4] = {1, 1, survivors, 1};
	int is_float = sizeof(Real)==4;

	int fpzip_bits = 8*sizeof(Real);
	if(getenv("FPZIP_BITS")) fpzip_bits = atoi(getenv("FPZIP_BITS"));
	if (fpzip_bits > 8*sizeof(Real)) fpzip_bits = 8*sizeof(Real);

//	fpz_decompress3D((char *)bufcompression + BITSETSIZE, survbytes, layout, tmp, (unsigned int *)&outbytes, 1, 8*sizeof(Real));
	fpz_decompress3D((char *)bufcompression + BITSETSIZE, survbytes, layout, tmp, (unsigned int *)&outbytes, is_float, fpzip_bits);
#if VERBOSE
	printf("fpzip3(D) : %d -> %d\n", survbytes, outbytes);
#endif

#elif defined(_USE_ZFP3_)
	double zfp_acc = 0.01;
	if(getenv("ZFP_ACC")) zfp_acc = atof(getenv("ZFP_ACC"));
	int is_float = sizeof(Real)==4;
	int layout[4] = {survivors, 1, 1, 1};
	size_t nbytes_zfp;
	int status = zfp_decompress_buffer((char *)((char *)bufcompression + BITSETSIZE), layout[0], layout[1], layout[2], zfp_acc, is_float, (unsigned char *)tmp, survbytes, &nbytes_zfp);
	outbytes = nbytes_zfp;
#if VERBOSE
	printf("zfp3(D) : %d -> %ld\n", survbytes, nbytes_zfp);
#endif
//	if (nbytes_zfp > survbytes) exit(1);

#elif defined(_USE_SZ3_)
	int layout[5] = {survivors, 1, 1, 0, 0};
	//int SZ_decompress_args(int dataType, unsigned char *bytes, int byteLength, void* decompressed_array, int r5, int r4, int r3, int r2, int r1);
        int sz_decompressedbytes = SZ_decompress_args(SZ_FLOAT, (unsigned char *)((char *)bufcompression + BITSETSIZE), survbytes, tmp, layout[4], layout[3], layout[2], layout[1], layout[0]);
	outbytes = sz_decompressedbytes;
//	outbytes = survivors*sizeof(float); 

	printf("sz3(D) : %d -> %d (vs %d)\n", survbytes, outbytes, survivors*sizeof(Real));
#endif

	memcpy((char *)bufcompression + BITSETSIZE, tmp, outbytes);
	free(tmp);

	bytes = BITSETSIZE + outbytes;

	}

#elif defined(_USE_SHUFFLE3_)

 #if !defined(_USE_FLOAT16_)	// SHUFFLE
	reshuffle3((char *)(bufcompression + BITSETSIZE), bytes - BITSETSIZE, sizeof(DataType));

 #else // _USE_FLOAT16_

	{
	reshuffle3((char *)(bufcompression + BITSETSIZE), bytes - BITSETSIZE, sizeof(unsigned short));

	size_t bytes_read = BITSETSIZE;
	
	//const
	int nelements = (bytes - bytes_read) / sizeof(unsigned short);
	if (nelements - expected == 1) {
		nelements--;
	}
	assert(expected == nelements);
	
	vector<DataType> datastream(nelements);
	
	{
		/* this is buggy if datastream is not 2-bytes aligned, not? */
		unsigned short elements[nelements];
		memcpy((void *)&elements, bufcompression + bytes_read, sizeof(unsigned short) * nelements);	
		for(int i = 0; i < nelements; ++i)
			datastream[i] =  _cvtfromf16(elements[i]);
	}
	
	full.load(datastream, mask);
	full.iwt(wtype);
	return;
	}
 #endif

#endif

	assert((bytes - sizeof(bitset<BS3>)) % sizeof(DataType) == 0 || float16);
	assert((bytes - sizeof(bitset<BS3>)) % sizeof(unsigned short) == 0);
	
//	bitset<BS3> mask;
//	const int expected = deserialize_bitset<BS3>(mask, bufcompression, BITSETSIZE);
	
	size_t bytes_read = BITSETSIZE;
	
	//const
	int nelements = (bytes - bytes_read) / (float16 ? sizeof(unsigned short) : sizeof(DataType));
	if (nelements - expected == 1) {
		nelements--;
	}
	assert(expected == nelements);
	
	vector<DataType> datastream(nelements);
	
	if (!float16)
	{
		memcpy((void *)&datastream.front(), bufcompression + bytes_read, sizeof(DataType) * nelements);	
	}
	else
	{
		/* this is buggy if datastream is not 2-bytes aligned, not? */
		unsigned short elements[nelements];
		memcpy((void *)&elements, bufcompression + bytes_read, sizeof(unsigned short) * nelements);	
		for(int i = 0; i < nelements; ++i)
			datastream[i] =  _cvtfromf16(elements[i]);
	}
	
	full.load(datastream, mask);
	full.iwt(wtype);
}

#ifdef _BLOCKSIZE_
template class WaveletCompressorGeneric<_BLOCKSIZE_, Real>;
template class WaveletCompressorGeneric_zlib<_BLOCKSIZE_, Real>;
#endif

#ifdef _VOXELS_ //mammamia whattahack
template class WaveletCompressorGeneric<_VOXELS_, float>;
template class WaveletCompressorGeneric_zlib<_VOXELS_, float>;
#endif
