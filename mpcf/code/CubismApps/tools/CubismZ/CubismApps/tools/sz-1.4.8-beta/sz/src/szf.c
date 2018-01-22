/**
 *  @file szf.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief the key C binding file to connect Fortran and C
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "sz.h"


//special notice: all the function names in this file must be lower-cases!!
void sz_init_c_(char *configFile,int *len,int *ierr)
{
    int i;
    char s2[*len+1];
    for(i=0;i<*len;i++)
        s2[i]=configFile[i];
    s2[*len]='\0';
 //   printf("sconfigFile=%s\n",configFile);
    *ierr = SZ_Init(s2);
}

void sz_finalize_c_()
{
	SZ_Finalize();
}

//compress with config (without args in function)
void sz_compress_d1_float_(float* data, unsigned char *bytes, int *outSize, int *r1)	
{
	unsigned char *tmp_bytes = SZ_compress(SZ_FLOAT, data, outSize, 0, 0, 0, 0, *r1);
	memcpy(bytes, tmp_bytes, *outSize);	
	free(tmp_bytes);
}

void sz_compress_d1_float_rev_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *r1)	
{
	unsigned char *tmp_bytes = SZ_compress_rev(SZ_FLOAT, data, reservedValue, outSize, 0, 0, 0, 0, *r1);
	memcpy(bytes, tmp_bytes, *outSize);	
	free(tmp_bytes);
}

void sz_compress_d2_float_(float* data, unsigned char *bytes, int *outSize, int *r1, int *r2)
{
	unsigned char *tmp_bytes = SZ_compress(SZ_FLOAT, data, outSize, 0, 0, 0, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d2_float_rev_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2)
{
	unsigned char *tmp_bytes = SZ_compress_rev(SZ_FLOAT, data, reservedValue, outSize, 0, 0, 0, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d3_float_(float* data, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3)
{
	unsigned char *tmp_bytes = SZ_compress(SZ_FLOAT, data, outSize, 0, 0, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d3_float_rev_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3)
{
	unsigned char *tmp_bytes = SZ_compress_rev(SZ_FLOAT, data, reservedValue, outSize, 0, 0, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d4_float_(float* data, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4)
{
	unsigned char *tmp_bytes = SZ_compress(SZ_FLOAT, data, outSize, 0, *r4, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d4_float_rev_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4)
{
	unsigned char *tmp_bytes = SZ_compress_rev(SZ_FLOAT, data, reservedValue, outSize, 0, *r4, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d5_float_(float* data, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4, int *r5)
{
	unsigned char *tmp_bytes = SZ_compress(SZ_FLOAT, data, outSize, *r5, *r4, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d5_float_rev_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4, int *r5)
{
	unsigned char *tmp_bytes = SZ_compress_rev(SZ_FLOAT, data, reservedValue, outSize, *r5, *r4, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d1_double_(double* data, unsigned char *bytes, int *outSize, int *r1)
{
	unsigned char *tmp_bytes = SZ_compress(SZ_DOUBLE, data, outSize, 0, 0, 0, 0, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d1_double_rev_(double* data, double *reservedValue, unsigned char *bytes, int *outSize, int *r1)
{
	unsigned char *tmp_bytes = SZ_compress_rev(SZ_DOUBLE, data, reservedValue, outSize, 0, 0, 0, 0, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d2_double_(double* data, unsigned char *bytes, int *outSize, int *r1, int *r2)
{
	unsigned char *tmp_bytes = SZ_compress(SZ_DOUBLE, data, outSize, 0, 0, 0, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d2_double_rev_(double* data, double *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2)
{
	unsigned char *tmp_bytes = SZ_compress_rev(SZ_DOUBLE, data, reservedValue, outSize, 0, 0, 0, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d3_double_(double* data, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3)
{
	unsigned char *tmp_bytes = SZ_compress(SZ_DOUBLE, data, outSize, 0, 0, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d3_double_rev_(double* data, double *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3)
{
	unsigned char *tmp_bytes = SZ_compress_rev(SZ_DOUBLE, data, reservedValue, outSize, 0, 0, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d4_double_(double* data, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4)
{
	unsigned char *tmp_bytes = SZ_compress(SZ_DOUBLE, data, outSize, 0, *r4, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d4_double_rev_(double* data, double *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4)
{
	unsigned char *tmp_bytes = SZ_compress_rev(SZ_DOUBLE, data, reservedValue, outSize, 0, *r4, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d5_double_(double* data, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4, int *r5)
{
	unsigned char *tmp_bytes = SZ_compress(SZ_DOUBLE, data, outSize, *r5, *r4, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d5_double_rev_(double* data, double *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4, int *r5)
{
	unsigned char *tmp_bytes = SZ_compress_rev(SZ_DOUBLE, data, reservedValue, outSize, *r5, *r4, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

//compress with args

void sz_compress_d1_float_args_(float* data, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1)
{
	unsigned char *tmp_bytes = SZ_compress_args(SZ_FLOAT, data, outSize, *errBoundMode, *absErrBound, *relBoundRatio, 0, 0, 0, 0, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d2_float_args_(float* data, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2)
{
	unsigned char *tmp_bytes = SZ_compress_args(SZ_FLOAT, data, outSize, *errBoundMode, *absErrBound, *relBoundRatio, 0, 0, 0, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d3_float_args_(float* data, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3)
{
	unsigned char *tmp_bytes = SZ_compress_args(SZ_FLOAT, data, outSize, *errBoundMode, *absErrBound, *relBoundRatio, 0, 0, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d4_float_args_(float* data, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4)
{
	unsigned char *tmp_bytes = SZ_compress_args(SZ_FLOAT, data, outSize, *errBoundMode, *absErrBound, *relBoundRatio, 0, *r4, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d5_float_args_(float* data, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5)
{
	unsigned char *tmp_bytes = SZ_compress_args(SZ_FLOAT, data, outSize, *errBoundMode, *absErrBound, *relBoundRatio, *r5, *r4, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d1_double_args_(double* data, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1)
{
	unsigned char *tmp_bytes = SZ_compress_args(SZ_DOUBLE, data, outSize, *errBoundMode, *absErrBound, *relBoundRatio, 0, 0, 0, 0, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d2_double_args_(double* data, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2)
{
	unsigned char *tmp_bytes = SZ_compress_args(SZ_DOUBLE, data, outSize, *errBoundMode, *absErrBound, *relBoundRatio, 0, 0, 0, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d3_double_args_(double* data, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3)
{
	unsigned char *tmp_bytes = SZ_compress_args(SZ_DOUBLE, data, outSize, *errBoundMode, *absErrBound, *relBoundRatio, 0, 0, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d4_double_args_(double* data, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4)
{
	unsigned char *tmp_bytes = SZ_compress_args(SZ_DOUBLE, data, outSize, *errBoundMode, *absErrBound, *relBoundRatio, 0, *r4, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d5_double_args_(double* data, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5)
{
	unsigned char *tmp_bytes = SZ_compress_args(SZ_DOUBLE, data, outSize, *errBoundMode, *absErrBound, *relBoundRatio, *r5, *r4, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

//--------------

void sz_compress_d1_float_rev_args_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1)
{
	unsigned char *tmp_bytes = SZ_compress_rev_args(SZ_FLOAT, data, reservedValue, outSize, *errBoundMode, *absErrBound, *relBoundRatio, 0, 0, 0, 0, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d2_float_rev_args_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2)
{
	unsigned char *tmp_bytes = SZ_compress_rev_args(SZ_FLOAT, data, reservedValue, outSize, *errBoundMode, *absErrBound, *relBoundRatio, 0, 0, 0, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d3_float_rev_args_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3)
{
	unsigned char *tmp_bytes = SZ_compress_rev_args(SZ_FLOAT, data, reservedValue, outSize, *errBoundMode, *absErrBound, *relBoundRatio, 0, 0, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d4_float_rev_args_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4)
{
	unsigned char *tmp_bytes = SZ_compress_rev_args(SZ_FLOAT, data, reservedValue, outSize, *errBoundMode, *absErrBound, *relBoundRatio, 0, *r4, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d5_float_rev_args_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5)
{
	unsigned char *tmp_bytes = SZ_compress_rev_args(SZ_FLOAT, data, reservedValue, outSize, *errBoundMode, *absErrBound, *relBoundRatio, *r5, *r4, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d1_double_rev_args_(double* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1)
{
	unsigned char *tmp_bytes = SZ_compress_rev_args(SZ_DOUBLE, data, reservedValue, outSize, *errBoundMode, *absErrBound, *relBoundRatio, 0, 0, 0, 0, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d2_double_rev_args_(double* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2)
{
	unsigned char *tmp_bytes = SZ_compress_rev_args(SZ_DOUBLE, data, reservedValue, outSize, *errBoundMode, *absErrBound, *relBoundRatio, 0, 0, 0, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d3_double_rev_args_(double* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3)
{
	unsigned char *tmp_bytes = SZ_compress_rev_args(SZ_DOUBLE, data, reservedValue, outSize, *errBoundMode, *absErrBound, *relBoundRatio, 0, 0, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
}

void sz_compress_d4_double_rev_args_(double* data, double *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4)
{
	unsigned char *tmp_bytes = SZ_compress_rev_args(SZ_DOUBLE, data, reservedValue, outSize, *errBoundMode, *absErrBound, *relBoundRatio, 0, *r4, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

void sz_compress_d5_double_rev_args_(double* data, double *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5)
{
	unsigned char *tmp_bytes = SZ_compress_rev_args(SZ_DOUBLE, data, reservedValue, outSize, *errBoundMode, *absErrBound, *relBoundRatio, *r5, *r4, *r3, *r2, *r1);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}

//decompress

void sz_decompress_d1_float_(unsigned char *bytes, int *byteLength, float *data, int *r1)
{
	float *tmp_data = SZ_decompress(SZ_FLOAT, bytes, *byteLength, 0, 0, 0, 0, *r1);
	memcpy(data, tmp_data, (*r1)*sizeof(float));
	free(tmp_data);
}

void sz_decompress_d2_float_(unsigned char *bytes, int *byteLength, float *data, int *r1, int *r2)
{
	int r;
	float *tmp_data = SZ_decompress(SZ_FLOAT, bytes, *byteLength, 0, 0, 0, *r2, *r1);
	r=(*r1)*(*r2);
	memcpy(data, tmp_data, r*sizeof(float));
	free(tmp_data);
}

void sz_decompress_d3_float_(unsigned char *bytes, int *byteLength, float *data, int *r1, int *r2, int *r3)
{
	int r;
	float *tmp_data = SZ_decompress(SZ_FLOAT, bytes, *byteLength, 0, 0, *r3, *r2, *r1);
	r=(*r1)*(*r2)*(*r3);
	memcpy(data, tmp_data, r*sizeof(float));
	free(tmp_data);
}

void sz_decompress_d4_float_(unsigned char *bytes, int *byteLength, float *data, int *r1, int *r2, int *r3, int *r4)
{
	int r;
	float *tmp_data = SZ_decompress(SZ_FLOAT, bytes, *byteLength, 0, *r4, *r3, *r2, *r1);
	r=(*r1)*(*r2)*(*r3)*(*r4);
	memcpy(data, tmp_data, r*sizeof(float));
	free(tmp_data);
}

void sz_decompress_d5_float_(unsigned char *bytes, int *byteLength, float *data, int *r1, int *r2, int *r3, int *r4, int *r5)
{
	int r;
	float *tmp_data = SZ_decompress(SZ_FLOAT, bytes, *byteLength, *r5, *r4, *r3, *r2, *r1);
	r=(*r1)*(*r2)*(*r3)*(*r4)*(*r5);
	memcpy(data, tmp_data, r*sizeof(float));
	free(tmp_data);
}

void sz_decompress_d1_double_(unsigned char *bytes, int *byteLength, double *data, int *r1)
{
	double *tmp_data = SZ_decompress(SZ_DOUBLE, bytes, *byteLength, 0, 0, 0, 0, *r1);
	memcpy(data, tmp_data, (*r1)*sizeof(double));
	free(tmp_data);
}

void sz_decompress_d2_double_(unsigned char *bytes, int *byteLength, double *data, int *r1, int *r2)
{
	int r;
	double *tmp_data = SZ_decompress(SZ_DOUBLE, bytes, *byteLength, 0, 0, 0, *r2, *r1);
	r=(*r1)*(*r2);
	memcpy(data, tmp_data, r*sizeof(double));
	free(tmp_data);
}

void sz_decompress_d3_double_(unsigned char *bytes, int *byteLength, double *data, int *r1, int *r2, int *r3)
{
	int r;
	double *tmp_data = SZ_decompress(SZ_DOUBLE, bytes, *byteLength, 0, 0, *r3, *r2, *r1);
	r=(*r1)*(*r2)*(*r3);
	memcpy(data, tmp_data, r*sizeof(double));
	free(tmp_data);
}

void sz_decompress_d4_double_(unsigned char *bytes, int *byteLength, double *data, int *r1, int *r2, int *r3, int *r4)
{
	int r;
	double *tmp_data = SZ_decompress(SZ_DOUBLE, bytes, *byteLength, 0, *r4, *r3, *r2, *r1);
	r=(*r1)*(*r2)*(*r3)*(*r4);
	memcpy(data, tmp_data, r*sizeof(double));
	free(tmp_data);
}

void sz_decompress_d5_double_(unsigned char *bytes, int *byteLength, double *data, int *r1, int *r2, int *r3, int *r4, int *r5)
{
	int r;
	double *tmp_data = SZ_decompress(SZ_DOUBLE, bytes, *byteLength, *r5, *r4, *r3, *r2, *r1);
	r=(*r1)*(*r2)*(*r3)*(*r4)*(*r5);
	memcpy(data, tmp_data, r*sizeof(double));
	free(tmp_data);
}

//-----------------TODO: batch mode-----------
void sz_batchaddvar_d1_float_(char* varName, int *len, float* data, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1)
{
	int i;
    char s2[*len+1];
    for(i=0;i<*len;i++)
        s2[i]=varName[i];
    s2[*len]='\0';		
	SZ_batchAddVar(s2, SZ_FLOAT, data, *errBoundMode, *absErrBound, *relBoundRatio, 0, 0, 0, 0, *r1);
}
void sz_batchaddvar_d2_float_(char* varName, int *len, float* data, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2)
{
	int i;
    char s2[*len+1];
    for(i=0;i<*len;i++)
        s2[i]=varName[i];
    s2[*len]='\0';		
	SZ_batchAddVar(s2, SZ_FLOAT, data, *errBoundMode, *absErrBound, *relBoundRatio, 0, 0, 0, *r2, *r1);
}
void sz_batchaddvar_d3_float_(char* varName, int *len, float* data, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3)
{
	int i;
    char s2[*len+1];
    for(i=0;i<*len;i++)
        s2[i]=varName[i];
    s2[*len]='\0';		
	SZ_batchAddVar(s2, SZ_FLOAT, data, *errBoundMode, *absErrBound, *relBoundRatio, 0, 0, *r3, *r2, *r1);
}
void sz_batchaddvar_d4_float_(char* varName, int *len, float* data, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4)
{
	int i;
    char s2[*len+1];
    for(i=0;i<*len;i++)
        s2[i]=varName[i];
    s2[*len]='\0';		
	SZ_batchAddVar(s2, SZ_FLOAT, data, *errBoundMode, *absErrBound, *relBoundRatio, 0, *r4, *r3, *r2, *r1);
}
void sz_batchaddvar_d5_float_(char* varName, int *len, float* data, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5)
{
	int i;
    char s2[*len+1];
    for(i=0;i<*len;i++)
        s2[i]=varName[i];
    s2[*len]='\0';		
	SZ_batchAddVar(s2, SZ_FLOAT, data, *errBoundMode, *absErrBound, *relBoundRatio, *r5, *r4, *r3, *r2, *r1);
}
void sz_batchaddvar_d1_double_(char* varName, int *len, double* data, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1)
{
	int i;
    char s2[*len+1];
    for(i=0;i<*len;i++)
        s2[i]=varName[i];
    s2[*len]='\0';		
	SZ_batchAddVar(s2, SZ_DOUBLE, data, *errBoundMode, *absErrBound, *relBoundRatio, 0, 0, 0, 0, *r1);
}
void sz_batchaddvar_d2_double_(char* varName, int *len, double* data, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2)
{
	int i;
    char s2[*len+1];
    for(i=0;i<*len;i++)
        s2[i]=varName[i];
    s2[*len]='\0';		
	SZ_batchAddVar(s2, SZ_DOUBLE, data, *errBoundMode, *absErrBound, *relBoundRatio, 0, 0, 0, *r2, *r1);
}
void sz_batchaddvar_d3_double_(char* varName, int *len, double* data, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3)
{
	int i;
    char s2[*len+1];
    for(i=0;i<*len;i++)
        s2[i]=varName[i];
    s2[*len]='\0';		
	SZ_batchAddVar(s2, SZ_DOUBLE, data, *errBoundMode, *absErrBound, *relBoundRatio, 0, 0, *r3, *r2, *r1);
}
void sz_batchaddvar_d4_double_(char* varName, int *len, double* data, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4)
{
	int i;
    char s2[*len+1];
    for(i=0;i<*len;i++)
        s2[i]=varName[i];
    s2[*len]='\0';		
	SZ_batchAddVar(s2, SZ_DOUBLE, data, *errBoundMode, *absErrBound, *relBoundRatio, 0, *r4, *r3, *r2, *r1);
}
void sz_batchaddvar_d5_double_(char* varName, int *len, double* data, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5)
{
	int i;
    char s2[*len+1];
    for(i=0;i<*len;i++)
        s2[i]=varName[i];
    s2[*len]='\0';		
	SZ_batchAddVar(s2, SZ_DOUBLE, data, *errBoundMode, *absErrBound, *relBoundRatio, *r5, *r4, *r3, *r2, *r1);
}
void sz_batchdelvar_c_(char* varName, int *len, int *errState)
{
	int i;
    char s2[*len+1];
    for(i=0;i<*len;i++)
        s2[i]=varName[i];
    s2[*len]='\0';
	*errState = SZ_batchDelVar(s2);
}
void sz_batch_compress_c_(unsigned char* bytes, int *outSize)
{
	unsigned char* tmp_bytes = SZ_batch_compress(outSize);
	memcpy(bytes, tmp_bytes, *outSize);
	free(tmp_bytes);
}
void sz_batch_decompress_c_(unsigned char* bytes, int *byteLength)
{
	SZ_batch_decompress(bytes, *byteLength);
}
void sz_getvardata_float_(char* varName, int *len, float* data, int *r1, int *r2, int *r3, int *r4, int *r5)
{
	int i;
    char s2[*len+1];
    for(i=0;i<*len;i++)
        s2[i]=varName[i];
    s2[*len]='\0';	
	
	float* tmp_data = (float*)SZ_getVarData(s2, r5, r4, r3, r2, r1);
	int size = computeDataLength(*r5, *r4, *r3, *r2, *r1);
	memcpy(data, tmp_data, size*sizeof(float));
	free(tmp_data);	
}
void sz_getvardata_double_(char* varName, int *len, double* data, int *r1, int *r2, int *r3, int *r4, int *r5)
{
	int i;
    char s2[*len+1];
    for(i=0;i<*len;i++)
        s2[i]=varName[i];
    s2[*len]='\0';	
    
	double* tmp_data = (double*)SZ_getVarData(s2, r5, r4, r3, r2, r1);
	int size = computeDataLength(*r5, *r4, *r3, *r2, *r1);
	memcpy(data, tmp_data, size*sizeof(double));
	free(tmp_data);
}


