/**
 *  @file TightDataPointStorageF.h
 *  @author Sheng Di and Dingwen Tao
 *  @date Aug, 2016
 *  @brief Header file for the tight data point storage (TDPS).
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _TightDataPointStorageF_H
#define _TightDataPointStorageF_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct TightDataPointStorageF
{
	int dataSeriesLength;
	int allSameData;
	double realPrecision;
	float medianValue;
	char reqLength;
	
	int exactDataNum;
	float reservedValue;
	
	unsigned char* rtypeArray;
	int rtypeArray_size;
	
	unsigned char* typeArray; //its size is dataSeriesLength/4 (or xxx/4+1) 
	int typeArray_size;
	
	unsigned char* leadNumArray; //its size is exactDataNum/4 (or exactDataNum/4+1)
	int leadNumArray_size;
	
	unsigned char* exactMidBytes;
	int exactMidBytes_size;
	
	unsigned char* residualMidBits;
	int residualMidBits_size;
	
	unsigned int intervals;
	
	unsigned char isLossless; //a mark to denote whether it's lossless compression (1 is yes, 0 is no)
} TightDataPointStorageF;

void new_TightDataPointStorageF_Empty(TightDataPointStorageF **this);
void new_TightDataPointStorageF_fromFlatBytes(TightDataPointStorageF **this, unsigned char* flatBytes, int flatBytesLength);
void decompressDataSeries_float_1D(float** data, int dataSeriesLength, TightDataPointStorageF* tdps);
void decompressDataSeries_float_2D(float** data, int r1, int r2, TightDataPointStorageF* tdps);
void decompressDataSeries_float_3D(float** data, int r1, int r2, int r3, TightDataPointStorageF* tdps);
void getSnapshotData_float_1D(float** data, int dataSeriesLength, TightDataPointStorageF* tdps);
void getSnapshotData_float_2D(float** data, int r1, int r2, TightDataPointStorageF* tdps);
void getSnapshotData_float_3D(float** data, int r1, int r2, int r3, TightDataPointStorageF* tdps);
void new_TightDataPointStorageF(TightDataPointStorageF **this, 
		int dataSeriesLength, int exactDataNum, 
		int* type, unsigned char* exactMidBytes, int exactMidBytes_size,
		unsigned char* leadNumIntArray,  //leadNumIntArray contains readable numbers....
		unsigned char* resiMidBits, int resiMidBits_size,
		unsigned char* resiBitLength, int resiBitLengthSize, 
		double realPrecision, float medianValue, char reqLength, unsigned int intervals);
void convertTDPStoFlatBytes_float(TightDataPointStorageF *tdps, unsigned char** bytes, int *size);
void free_TightDataPointStorageF(TightDataPointStorageF *tdps);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _TightDataPointStorageF_H  ----- */
