#ifndef _MYFPZIP_H_
#define _MYFPZIP_H_ 1

#include "fpzip.h"
#include <climits>

#if defined(_FPZIP_24_)
#define FPZIP_ACC 24
#elif defined(_FPZIP_16_)
#define FPZIP_ACC 16
#else
#define FPZIP_ACC 32
#endif

static void fpz_compress3D(void *in, unsigned int inbytes, int layout[4], void *out, unsigned int *outbytes, int isfloat, int prec)
{
	int i;
	int precbits = prec; 
	int type = isfloat ? FPZIP_TYPE_FLOAT : FPZIP_TYPE_DOUBLE;
	int dim[4] = {32, 32, 32, 1};

	for (i = 0; i < 3; i++) dim[i] = layout[i];

	FPZ* fpz = fpzip_write_to_buffer(out, inbytes /*+1024*/);
	fpz->type = type;
	fpz->prec = precbits;
//	fpz->nx = dim[0];
//	fpz->ny = dim[1];
//	fpz->nz = dim[2];
//	fpz->nf = 1;

	fpz->nf = dim[0];
	fpz->nz = dim[1];
	fpz->ny = dim[2];
	fpz->nx = 1;

	// perform actual compression
	*outbytes = fpzip_write(fpz, in);
	if (*outbytes == 0) {
		fprintf(stderr, "compression failed: %s\n", fpzip_errstr[fpzip_errno]);
		return;
	}
	fpzip_write_close(fpz);

//	*outbytes = fpzip_memory_write(out, inbytes, in, prec, dp, dim[0], dim[1], dim[2], dim[3]);
}

static void fpz_decompress3D(char *in, unsigned int inbytes, int layout[4], char *out, unsigned int *outbytes, int isfloat, int prec)
{
	int i;
	int precbits = prec; 
	int type = isfloat ? FPZIP_TYPE_FLOAT : FPZIP_TYPE_DOUBLE;
	int dim[4] = {32, 32, 32, 1};
	int targetbytes;

	for (i = 0; i < 3; i++) dim[i] = layout[i];
	//printf("fpz_decompress3D\n");
	if (isfloat) {
		targetbytes = dim[0]*dim[1]*dim[2]*dim[3]*sizeof(float);
	}
	else {
		targetbytes = dim[0]*dim[1]*dim[2]*dim[3]*sizeof(double);
	}

	FPZ* fpz = fpzip_read_from_buffer(in);
	fpz->type = type;
	fpz->prec = precbits;
//	fpz->nx = dim[0];
//	fpz->ny = dim[1];
//	fpz->nz = dim[2];
//	fpz->nf = 1;

	fpz->nf = dim[0];
	fpz->nz = dim[1];
	fpz->ny = dim[2];
	fpz->nx = 1;

	// perform actual compression
	if (!fpzip_read(fpz, out)) {
		fprintf(stderr, "decompression failed: %s\n", fpzip_errstr[fpzip_errno]);
		*outbytes = -1;
		return;
	}

//	readbytes = fpzip_memory_read(in, out, prec, dp, dim[0], dim[1], dim[2], dim[3]);
	
	*outbytes = targetbytes;
	return;
}

#endif
