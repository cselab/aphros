/**
 *  @file TightDataPointStorageD.h
 *  @author Sheng Di
 *  @date April, 2016
 *  @brief Header file for the tight data point storage (TDPS).
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _TightDataPointStorageD_H
#define _TightDataPointStorageD_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct TightDataPointStorageD
{
	int dataSeriesLength;
	int allSameData;
	int bestChoice;
	double realPrecision;
	double medianValue;
	char reqLength;	
	
	int exactDataNum;
	double reservedValue;
	
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
	
	unsigned char isLossless;

} TightDataPointStorageD;

void new_TightDataPointStorageD_Empty(TightDataPointStorageD **this);
void new_TightDataPointStorageD_fromFlatBytes(TightDataPointStorageD **this, unsigned char* flatBytes, int flatBytesLength);
void decompressDataSeries_double_1D(double** data, int dataSeriesLength, TightDataPointStorageD* tdps);
void decompressDataSeries_double_2D(double** data, int r1, int r2, TightDataPointStorageD* tdps);
void decompressDataSeries_double_3D(double** data, int r1, int r2, int r3, TightDataPointStorageD* tdps);
void getSnapshotData_double_1D(double** data, int dataSeriesLength, TightDataPointStorageD* tdps);
void getSnapshotData_double_2D(double** data, int r1, int r2, TightDataPointStorageD* tdps);
void getSnapshotData_double_3D(double** data, int r1, int r2, int r3, TightDataPointStorageD* tdps);
void new_TightDataPointStorageD(TightDataPointStorageD **this, 
		int dataSeriesLength, int exactDataNum, 
		int* type, unsigned char* exactMidBytes, int exactMidBytes_size,
		unsigned char* leadNumIntArray,  //leadNumIntArray contains readable numbers....
		unsigned char* resiMidBits, int resiMidBits_size,
		unsigned char* resiBitLength, int resiBitLengthSize, 
		double realPrecision, double medianValue, char reqLength, unsigned int intervals);
void convertTDPStoFlatBytes_double(TightDataPointStorageD *tdps, unsigned char** bytes, int *size);
void free_TightDataPointStorageD(TightDataPointStorageD *tdps);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _TightDataPointStorageD_H  ----- */
