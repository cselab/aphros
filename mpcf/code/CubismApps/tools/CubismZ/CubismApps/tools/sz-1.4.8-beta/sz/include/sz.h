/**
 *  @file sz.h
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief Header file for the whole detector.
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _SZ_H
#define _SZ_H

#include <stdio.h>
#include <sys/time.h>      /* For gettimeofday(), in microseconds */
#include <time.h>          /* For time(), in seconds */
#include "iniparser.h"
#include "CompressElement.h"
#include "DynamicByteArray.h"
#include "DynamicIntArray.h"
#include "VarSet.h"
#include "Huffman.h"

#ifdef _WIN32
#define PATH_SEPARATOR ';'
#else
#define PATH_SEPARATOR ':'
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define SZ_VERNUM 0x0130
#define SZ_VER_MAJOR 1
#define SZ_VER_MINOR 4
#define SZ_VER_REVISION 8

#define HZ 102
#define SZ 101

#define ABS 0
#define REL 1
#define ABS_AND_REL 2
#define ABS_OR_REL 3

#define SZ_FLOAT 0
#define SZ_DOUBLE 1

#define LITTLE_ENDIAN_DATA 0
#define BIG_ENDIAN_DATA 1 //big_endian (ppc, max, etc.) ; little_endian (x86, x64, etc.)

#define LITTLE_ENDIAN_SYSTEM 0
#define BIG_ENDIAN_SYSTEM 1

#define DynArrayInitLen 1024

#define MIN_ZLIB_DEC_ALLOMEM_BYTES 1000000

//#define maxRangeRadius 32768
//#define maxRangeRadius 1048576//131072

#define SZ_BEST_SPEED 0
#define SZ_BEST_COMPRESSION 1
#define SZ_DEFAULT_COMPRESSION 2

//Note: the following setting should be consistent with stateNum in Huffman.h
//#define intvCapacity 65536
//#define intvRadius 32768
//#define intvCapacity 131072
//#define intvRadius 65536

extern unsigned int maxRangeRadius;

extern int intvCapacity;
extern int intvRadius;

extern int sysEndianType; //endian type of the system
extern int dataEndianType; //endian type of the data
//extern int maxSegmentNum;

extern char maxHeap[10];
 
extern long status;

extern int sol_ID;
extern int errorBoundMode; //ABS, REL, ABS_AND_REL, or ABS_OR_REL

extern int gzipMode; //four options: Z_NO_COMPRESSION, or Z_BEST_SPEED, Z_BEST_COMPRESSION, Z_DEFAULT_COMPRESSION

extern char *sz_cfgFile;

extern int offset;

extern double absErrBound;
extern double relBoundRatio;

extern int versionNumber[3];

extern int layers;
extern float predThreshold;
extern int sampleDistance;
extern char optQuantMode;

extern int szMode; //0 (best speed) or 1 (better compression with Gzip)

//extern int spaceFillingCurveTransform; //default is 0, or 1 set by sz.config
//extern int reOrgSize; //the granularity of the reganization of the original data

extern SZ_VarSet* sz_varset;

//typedef unsigned long unsigned long;
//typedef unsigned int uint;

typedef union lshort
{
	unsigned short svalue;
	unsigned char byte[2];
} lshort;

typedef union ldouble
{
    double value;
    unsigned long lvalue;
    unsigned char byte[8];
} ldouble;

typedef union lfloat
{
    float value;
    unsigned int ivalue;
    unsigned char byte[4];
} lfloat;

/* array meta data and compression parameters for SZ_Init_Params() */
typedef struct sz_params
{
	unsigned int max_quant_intervals;
	unsigned int quantization_intervals;
    int dataEndianType;
    int sysEndianType; //sysEndianType is actually set automatically.
    int sol_ID;
    int layers;
    int sampleDistance;
    float predThreshold;    
    int offset;
    int szMode;
    int gzipMode;
    int  errorBoundMode;
    double absErrBound;
    double relBoundRatio;
} sz_params;

extern sz_params *conf_params;

//conf.c
void updateQuantizationInfo(int quant_intervals);
void clearHuffmanMem();
int SZ_ReadConf();
int SZ_LoadConf();
int checkVersion(char* version);
unsigned int roundUpToPowerOf2(unsigned int base);

//dataCompression.c
void computeRangeSize_double(double* oriData, int size, double* valueRangeSize, double* medianValue);
void computeRangeSize_float(float* oriData, int size, float* valueRangeSize, float* medianValue);
double min_d(double a, double b);
double max_d(double a, double b);
float min_f(float a, float b);
float max_f(float a, float b);
double getRealPrecision_double(double valueRangeSize, int errBoundMode, double absErrBound, double relBoundRatio);
double getRealPrecision_float(float valueRangeSize, int errBoundMode, double absErrBound, double relBoundRatio);
void symTransform_8bytes(unsigned char data[8]);
void flush_to_bigEndian_8bytes(unsigned char data[8], int dataEndianType);
void symTransform_2bytes(unsigned char data[2]);
void symTransform_4bytes(unsigned char data[4]);
void flush_to_bigEndian_4bytes(unsigned char data[4]);
void bigEndian_to_OSEndian_double(unsigned char data[8]);
void bigEndian_to_OSEndian_float(unsigned char data[4]);
void compressSingleFloatValue(FloatValueCompressElement *vce, float tgtValue, float precision, float medianValue, 
		int reqLength, int reqBytesLength, int resiBitsLength);
void compressSingleDoubleValue(DoubleValueCompressElement *vce, double tgtValue, double precision, double medianValue, 
		int reqLength, int reqBytesLength, int resiBitsLength);
int compIdenticalLeadingBytesCount_double(unsigned char* preBytes, unsigned char* curBytes);
int compIdenticalLeadingBytesCount_float(unsigned char* preBytes, unsigned char* curBytes);
void addExactData(DynamicByteArray *exactMidByteArray, DynamicIntArray *exactLeadNumArray, 
		DynamicIntArray *resiBitArray, LossyCompressionElement *lce);

//ByteToolkit.c
int bytesToInt_bigEndian(unsigned char* bytes);
void intToBytes_bigEndian(unsigned char *b, unsigned int num);
long bytesToLong_bigEndian(unsigned char* b);
void longToBytes_bigEndian(unsigned char *b, unsigned long num);
long doubleToOSEndianLong(double value);
int floatToOSEndianInt(float value);
short getExponent_float(float value);
short getPrecisionReqLength_float(float precision);
short getExponent_double(double value);
short getPrecisionReqLength_double(double precision);
unsigned char numberOfLeadingZeros_Int(int i);
unsigned char numberOfLeadingZeros_Long(long i);
unsigned char getLeadingNumbers_Int(int v1, int v2);
unsigned char getLeadingNumbers_Long(long v1, long v2);
short bytesToShort(unsigned char* bytes);
int bytesToInt(unsigned char* bytes);
long bytesToLong(unsigned char* bytes);
float bytesToFloat(unsigned char* bytes);
void floatToBytes(unsigned char *b, float num);
double bytesToDouble(unsigned char* bytes);
void doubleToBytes(unsigned char *b, double num);
int extractBytes(unsigned char* byteArray, int k, int validLength);
int getMaskRightCode(int m);
int getLeftMovingCode(int kMod8);
int getRightMovingSteps(int kMod8, int resiBitLength);
int getRightMovingCode(int kMod8, int resiBitLength);
unsigned short* convertByteDataToShortArray(unsigned char* bytes, int byteLength);
void convertShortArrayToBytes(unsigned short* states, int stateLength, unsigned char* bytes);

//TypeManager.c
int convertIntArray2ByteArray_fast_2b(unsigned char* timeStepType, int timeStepTypeLength, unsigned char **result);
void convertByteArray2IntArray_fast_2b(int stepLength, unsigned char* byteArray, int byteArrayLength, unsigned char **intArray);
int convertIntArray2ByteArray_fast_3b(unsigned char* timeStepType, int timeStepTypeLength, unsigned char **result);
void convertByteArray2IntArray_fast_3b(int stepLength, unsigned char* byteArray, int byteArrayLength, unsigned char **intArray);
int getLeftMovingSteps(int k, unsigned char resiBitLength);
int convertIntArray2ByteArray_fast_dynamic(unsigned char* timeStepType, unsigned char* resiBitLength, int resiBitLengthLength, unsigned char **bytes);
int computeBitNumRequired(int dataLength);
void decompressBitArraybySimpleLZ77(int** result, unsigned char* bytes, int bytesLength, int totalLength, int validLength);

//test_zlib.c
unsigned long zlib_compress(unsigned char* data, unsigned long dataLength, unsigned char** compressBytes, int level);
unsigned long zlib_compress2(unsigned char* data, unsigned long dataLength, unsigned char** compressBytes, int level);
unsigned long zlib_uncompress(unsigned char* compressBytes, unsigned long cmpSize, unsigned char** oriData, unsigned long targetOriSize);
unsigned long zlib_uncompress2(unsigned char* compressBytes, unsigned long cmpSize, unsigned char** oriData, unsigned long targetOriSize);

//szf.c
void sz_init_c_(char *configFile,int *len,int *ierr);
void sz_finalize_c_();
void SZ_writeData_inBinary_d1_Float_(float* data, char *fileName, int *len);
void sz_compress_d1_float_(float* data, unsigned char *bytes, int *outSize, int *r1);
void sz_compress_d1_float_rev_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *r1);
void sz_compress_d2_float_(float* data, unsigned char *bytes, int *outSize, int *r1, int *r2);
void sz_compress_d2_float_rev_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2);
void sz_compress_d3_float_(float* data, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3);
void sz_compress_d3_float_rev_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3);
void sz_compress_d4_float_(float* data, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d4_float_rev_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d5_float_(float* data, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4, int *r5);
void sz_compress_d5_float_rev_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4, int *r5);

void sz_compress_d1_double_(double* data, unsigned char *bytes, int *outSize, int *r1);
void sz_compress_d1_double_rev_(double* data, double *reservedValue, unsigned char *bytes, int *outSize, int *r1);
void sz_compress_d2_double_(double* data, unsigned char *bytes, int *outSize, int *r1, int *r2);
void sz_compress_d2_double_rev_(double* data, double *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2);
void sz_compress_d3_double_(double* data, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3);
void sz_compress_d3_double_rev_(double* data, double *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3);
void sz_compress_d4_double_(double* data, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d4_double_rev_(double* data, double *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d5_double_(double* data, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4, int *r5);
void sz_compress_d5_double_rev_(double* data, double *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4, int *r5);

void sz_compress_d1_float_args_(float* data, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1);
void sz_compress_d2_float_args_(float* data, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2);
void sz_compress_d3_float_args_(float* data, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3);
void sz_compress_d4_float_args_(float* data, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d5_float_args_(float* data, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5);
void sz_compress_d1_double_args_(double* data, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1);
void sz_compress_d2_double_args_(double* data, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2);
void sz_compress_d3_double_args_(double* data, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3);
void sz_compress_d4_double_args_(double* data, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d5_double_args_(double* data, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5);

void sz_compress_d1_float_rev_args_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1);
void sz_compress_d2_float_rev_args_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2);
void sz_compress_d3_float_rev_args_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3);
void sz_compress_d4_float_rev_args_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d5_float_rev_args_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5);
void sz_compress_d1_double_rev_args_(double* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1);
void sz_compress_d2_double_rev_args_(double* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2);
void sz_compress_d3_double_rev_args_(double* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3);
void sz_compress_d4_double_rev_args_(double* data, double *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d5_double_rev_args_(double* data, double *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5);

void sz_decompress_d1_float_(unsigned char *bytes, int *byteLength, float *data, int *r1);
void sz_decompress_d2_float_(unsigned char *bytes, int *byteLength, float *data, int *r1, int *r2);
void sz_decompress_d3_float_(unsigned char *bytes, int *byteLength, float *data, int *r1, int *r2, int *r3);
void sz_decompress_d4_float_(unsigned char *bytes, int *byteLength, float *data, int *r1, int *r2, int *r3, int *r4);
void sz_decompress_d5_float_(unsigned char *bytes, int *byteLength, float *data, int *r1, int *r2, int *r3, int *r4, int *r5);
void sz_decompress_d1_double_(unsigned char *bytes, int *byteLength, double *data, int *r1);
void sz_decompress_d2_double_(unsigned char *bytes, int *byteLength, double *data, int *r1, int *r2);
void sz_decompress_d3_double_(unsigned char *bytes, int *byteLength, double *data, int *r1, int *r2, int *r3);
void sz_decompress_d4_double_(unsigned char *bytes, int *byteLength, double *data, int *r1, int *r2, int *r3, int *r4);
void sz_decompress_d5_double_(unsigned char *bytes, int *byteLength, double *data, int *r1, int *r2, int *r3, int *r4, int *r5);

void sz_batchaddVar_d1_float_(char* varName, int *len, float* data, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1);
void sz_batchaddvar_d2_float_(char* varName, int *len, float* data, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2);
void sz_batchaddvar_d3_float_(char* varName, int *len, float* data, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3);
void sz_batchaddvar_d4_float_(char* varName, int *len, float* data, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4);
void sz_batchaddvar_d5_float_(char* varName, int *len, float* data, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5);
void sz_batchaddvar_d1_double_(char* varName, int *len, double* data, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1);
void sz_batchaddvar_d2_double_(char* varName, int *len, double* data, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2);
void sz_batchaddvar_d3_double_(char* varName, int *len, double* data, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3);
void sz_batchaddvar_d4_double_(char* varName, int *len, double* data, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4);
void sz_batchaddvar_d5_double_(char* varName, int *len, double* data, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5);
void sz_batchdelvar_c_(char* varName, int *len, int *errState);
void sz_batch_compress_c_(unsigned char* bytes, int *outSize);
void sz_batch_decompress_c_(unsigned char* bytes, int *byteLength);
void sz_getvardata_float_(char* varName, int *len, float* data, int *r1, int *r2, int *r3, int *r4, int *r5);
void sz_getvardata_double_(char* varName, int *len, double* data, int *r1, int *r2, int *r3, int *r4, int *r5);


//sz.h
void SZ_Reset();
int SZ_Init(char *configFilePath);
int SZ_Init_Params(sz_params *params);
int computeDataLength(int r5, int r4, int r3, int r2, int r1);
int computeDimension(int r5, int r4, int r3, int r2, int r1);

int getPredictionCoefficients(int layers, int dimension, int **coeff_array);
unsigned int optimize_intervals_float_1D(float *oriData, int dataLength, double realPrecision);
unsigned int optimize_intervals_float_2D(float *oriData, int r1, int r2, double realPrecision);
unsigned int optimize_intervals_float_3D(float *oriData, int r1, int r2, int r3, double realPrecision);
unsigned int optimize_intervals_double_1D(double *oriData, int dataLength, double realPrecision);
unsigned int optimize_intervals_double_2D(double *oriData, int r1, int r2, double realPrecision);
unsigned int optimize_intervals_double_3D(double *oriData, int r1, int r2, int r3, double realPrecision);

//void SZ_compress_args_float_NoCkRngeNoGzip(char** newByteData, float *oriData, int dataLength, double realPrecision, int *outSize);
void SZ_compress_args_float_NoCkRngeNoGzip_1D(unsigned char** newByteData, float *oriData, int dataLength, double realPrecision, int *outSize, float valueRangeSize, float medianValue_f);
void SZ_compress_args_float_NoCkRngeNoGzip_2D(unsigned char** newByteData, float *oriData, int r1, int r2, double realPrecision, int *outSize, float valueRangeSize, float medianValue_f);
void SZ_compress_args_float_NoCkRngeNoGzip_3D(unsigned char** newByteData, float *oriData, int r1, int r2, int r3, double realPrecision, int *outSize, float valueRangeSize, float medianValue_f);
void SZ_compress_args_double_NoCkRngeNoGzip_1D(unsigned char** newByteData, double *oriData, int dataLength, double realPrecision, int *outSize, double valueRangeSize, double medianValue_d);
void SZ_compress_args_double_NoCkRngeNoGzip_2D(unsigned char** newByteData, double *oriData, int r1, int r2, double realPrecision, int *outSize, double valueRangeSize, double medianValue_d);
void SZ_compress_args_double_NoCkRngeNoGzip_3D(unsigned char** newByteData, double *oriData, int r1, int r2, int r3, double realPrecision, int *outSize, double valueRangeSize, double medianValue_d);

void SZ_compress_args_float_withinRange(unsigned char** newByteData, float *oriData, int dataLength, int *outSize);
void SZ_compress_args_double_withinRange(unsigned char** newByteData, double *oriData, int dataLength, int *outSize);

void SZ_compress_args_float(unsigned char** newByteData, float *oriData, 
int r5, int r4, int r3, int r2, int r1, int *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio);
void SZ_compress_args_double(unsigned char** newByteData, double *oriData, 
int r5, int r4, int r3, int r2, int r1, int *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio);

void SZ_compress_args_float_wRngeNoGzip(unsigned char** newByteData, float *oriData, 
int r5, int r4, int r3, int r2, int r1, int *outSize, 
int errBoundMode, double absErr_Bound, double rel_BoundRatio);
void SZ_compress_args_double_wRngeNoGzip(unsigned char** newByteData, double *oriData, 
int r5, int r4, int r3, int r2, int r1, int *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio);

unsigned char *SZ_compress(int dataType, void *data, int *outSize, int r5, int r4, int r3, int r2, int r1);
unsigned char *SZ_compress_args(int dataType, void *data, int *outSize, int errBoundMode, double absErrBound, double relBoundRatio, int r5, int r4, int r3, int r2, int r1);
int SZ_compress_args2(int dataType, void *data, unsigned char* compressed_bytes, int *outSize, int errBoundMode, double absErrBound, double relBoundRatio, int r5, int r4, int r3, int r2, int r1);

unsigned char *SZ_compress_rev_args(int dataType, void *data, void *reservedValue, int *outSize, int errBoundMode, double absErrBound, double relBoundRatio, int r5, int r4, int r3, int r2, int r1);
int SZ_compress_rev_args2(int dataType, void *data, void *reservedValue, unsigned char* compressed_bytes, int *outSize, int errBoundMode, double absErrBound, double relBoundRatio, int r5, int r4, int r3, int r2, int r1);
unsigned char *SZ_compress_rev(int dataType, void *data, void *reservedValue, int *outSize, int r5, int r4, int r3, int r2, int r1);

void SZ_decompress_args_float(float** newData, int r5, int r4, int r3, int r2, int r1, unsigned char* cmpBytes, int cmpSize);
void SZ_decompress_args_double(double** newData, int r5, int r4, int r3, int r2, int r1, unsigned char* cmpBytes, int cmpSize);
void *SZ_decompress(int dataType, unsigned char *bytes, int byteLength, int r5, int r4, int r3, int r2, int r1);
int SZ_decompress_args(int dataType, unsigned char *bytes, int byteLength, void* decompressed_array, int r5, int r4, int r3, int r2, int r1);

void filloutDimArray(int* dim, int r5, int r4, int r3, int r2, int r1);
unsigned char* SZ_batch_compress(int *outSize);
SZ_VarSet* SZ_batch_decompress(unsigned char* compressedStream, int length);

void SZ_Finalize();

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _SZ_H  ----- */
