/**
 *  @file TightPointDataStorageD.c
 *  @author Sheng Di
 *  @date May, 2016
 *  @brief The functions used to construct the tightPointDataStorage element for storing compressed bytes.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include "TightDataPointStorageD.h"
#include "sz.h"
#include "Huffman.h"
//#include "rw.h"

void new_TightDataPointStorageD_Empty(TightDataPointStorageD **this)
{
	*this = (TightDataPointStorageD*)malloc(sizeof(TightDataPointStorageD));
	(*this)->dataSeriesLength = 0;
	(*this)->allSameData = 0;
	(*this)->exactDataNum = 0;
	(*this)->reservedValue = 0;
	
	(*this)->rtypeArray = NULL;
	(*this)->rtypeArray_size = 0;
	
	(*this)->typeArray = NULL; //its size is dataSeriesLength/4 (or xxx/4+1) 
	(*this)->typeArray_size = 0;
	
	(*this)->leadNumArray = NULL; //its size is exactDataNum/4 (or exactDataNum/4+1)
	(*this)->leadNumArray_size = 0;
	
	(*this)->exactMidBytes = NULL;
	(*this)->exactMidBytes_size = 0;
	
	(*this)->residualMidBits = NULL;
	(*this)->residualMidBits_size = 0;
	
	(*this)->intervals = 0;
}

void new_TightDataPointStorageD_fromFlatBytes(TightDataPointStorageD **this, unsigned char* flatBytes, int flatBytesLength)
{
	new_TightDataPointStorageD_Empty(this);
	int i, index = 0;
	char version[3];
	for (i = 0; i < 3; i++)
		version[i] = flatBytes[index++]; //3
	if(checkVersion(version)!=1)
	{
		//wrong version
		printf("Wrong version: \nCompressed-data version (%d.%d.%d)\n",version[0], version[1], version[2]);
		printf("Current sz version: (%d.%d.%d)\n", versionNumber[0], versionNumber[1], versionNumber[2]);
		exit(0);
	}
	
	unsigned char dsLengthBytes[4];
	for (i = 0; i < 4; i++)
		dsLengthBytes[i] = flatBytes[index++];
		
	//TODO
	(*this)->dataSeriesLength = bytesToInt_bigEndian(dsLengthBytes);// 4
	unsigned char sameRByte = flatBytes[index++]; //1
	int same = sameRByte & 0x01;
	szMode = (sameRByte & 0x06)>>1;
	(*this)->isLossless = (sameRByte & 0x16)>>4;
	//printf("szMode=%d\n",szMode);

	if((*this)->isLossless==1)
	{
		(*this)->exactMidBytes = flatBytes+8;
		return;
	}
	else if(same==1)
	{
		(*this)->allSameData = 1;
		int exactMidBytesLength = flatBytesLength - 3 - 4 -1;
		if(exactMidBytesLength>0)
			(*this)->exactMidBytes = (unsigned char*)malloc(sizeof(unsigned char)*exactMidBytesLength);
		else
			(*this)->exactMidBytes = NULL;
		for(i = 0;i<exactMidBytesLength;i++)
			(*this)->exactMidBytes[i] = flatBytes[index++];
		return;
	}
	else
		(*this)->allSameData = 0;
		
	int rtype_ = sameRByte & 0x08; //1000		

	unsigned char byteBuf[8];
	unsigned char byteBuf2[4];	
	
	for (i = 0; i < 4; i++)
		byteBuf2[i] = flatBytes[index++];
	(*this)->intervals = bytesToInt_bigEndian(byteBuf2);// 4	

	for (i = 0; i < 8; i++)
		byteBuf[i] = flatBytes[index++];
	(*this)->medianValue = bytesToDouble(byteBuf);//8

	(*this)->reqLength = flatBytes[index++]; //1

	for (i = 0; i < 8; i++)
		byteBuf[i] = flatBytes[index++];
	(*this)->realPrecision = bytesToDouble(byteBuf);//8

	for (i = 0; i < 4; i++)
		byteBuf2[i] = flatBytes[index++];
	(*this)->typeArray_size = bytesToInt_bigEndian(byteBuf2);// 4		

	if(rtype_!=0)
	{
		for(i = 0;i<4;i++) 
			byteBuf[i] = flatBytes[index++];
		(*this)->rtypeArray_size = bytesToInt_bigEndian(byteBuf);//4		
	}
	else
		(*this)->rtypeArray_size = 0;

	for (i = 0; i < 4; i++)
		byteBuf2[i] = flatBytes[index++];
	(*this)->exactDataNum = bytesToInt_bigEndian(byteBuf2);// 4

	for (i = 0; i < 4; i++)
		byteBuf2[i] = flatBytes[index++];
	(*this)->exactMidBytes_size = bytesToInt_bigEndian(byteBuf2);// 4

	int typeArrayLength = 0;
	if (rtype_ != 0) {
		if((*this)->rtypeArray_size>0)
			(*this)->rtypeArray = (unsigned char*)malloc(sizeof(unsigned char)*(*this)->rtypeArray_size);
		else
			(*this)->rtypeArray = NULL;
			
		for (i = 0; i < 8; i++)
			byteBuf[i] = flatBytes[index++];
		(*this)->reservedValue = bytesToDouble(byteBuf);//8
		
	}

	(*this)->typeArray = (unsigned char*)malloc(sizeof(unsigned char)*(*this)->typeArray_size);
		
	int logicLeadNumBitsNum = (*this)->exactDataNum * 2;
	if (logicLeadNumBitsNum % 8 == 0)
	{
		(*this)->leadNumArray_size = logicLeadNumBitsNum >> 3;
	}
	else
	{
		(*this)->leadNumArray_size = (logicLeadNumBitsNum >> 3) + 1;
	}
	if((*this)->leadNumArray_size>0)
		(*this)->leadNumArray = (unsigned char*)malloc(sizeof(unsigned char)*(*this)->leadNumArray_size);
	else
		(*this)->leadNumArray = NULL;
	
	if((*this)->exactMidBytes_size>0)
		(*this)->exactMidBytes = (unsigned char*)malloc(sizeof(unsigned char)*(*this)->exactMidBytes_size);
	else
		(*this)->exactMidBytes = NULL;

	if ((*this)->rtypeArray != NULL) 
	{
		(*this)->residualMidBits_size = flatBytesLength - 3 - 4 - 1 - 4 - 8 - 1 - 8 - 4 - 4 - 4 - 4
				- 8 - (*this)->rtypeArray_size 
				- (*this)->typeArray_size - (*this)->leadNumArray_size
				- (*this)->exactMidBytes_size;
		for (i = 0; i < (*this)->rtypeArray_size; i++)
			(*this)->rtypeArray[i] = flatBytes[index++];
	}
	else
	{
		(*this)->residualMidBits_size = flatBytesLength - 3 - 4 - 1 - 4 - 8 - 1 - 8 - 4 - 4
				- 4 - (*this)->typeArray_size
				- (*this)->leadNumArray_size - (*this)->exactMidBytes_size;
	}	
	if((*this)->residualMidBits_size>0)
		(*this)->residualMidBits = (unsigned char*)malloc(sizeof(unsigned char)*(*this)->residualMidBits_size);
	else
		(*this)->residualMidBits = NULL;
	//for (i = 0; i < (*this)->typeArray_size; i++)
	//	(*this)->typeArray[i] = flatBytes[index++];
	memcpy((*this)->typeArray, &flatBytes[index], (*this)->typeArray_size*sizeof(char));
	index+=(*this)->typeArray_size*sizeof(char);
	for (i = 0; i < (*this)->leadNumArray_size; i++)
		(*this)->leadNumArray[i] = flatBytes[index++];
	for (i = 0; i < (*this)->exactMidBytes_size; i++)
	{
		(*this)->exactMidBytes[i] = flatBytes[index++];
	}	
	for (i = 0; i < (*this)->residualMidBits_size; i++)
		(*this)->residualMidBits[i] = flatBytes[index++];	
}

void decompressDataSeries_double_1D(double** data, int dataSeriesLength, TightDataPointStorageD* tdps) 
{
	int i, j, k = 0, p = 0, l = 0; // k is to track the location of residual_bit
								// in resiMidBits, p is to track the
								// byte_index of resiMidBits, l is for
								// leadNum
	updateQuantizationInfo(tdps->intervals);

	unsigned char* leadNum;
	double interval = tdps->realPrecision*2;
	
	convertByteArray2IntArray_fast_2b(tdps->exactDataNum, tdps->leadNumArray, tdps->leadNumArray_size, &leadNum);
	*data = (double*)malloc(sizeof(double)*dataSeriesLength);

	int* type = (int*)malloc(dataSeriesLength*sizeof(int));
	//convertByteArray2IntArray_fast_3b(dataSeriesLength, tdps->typeArray, tdps->typeArray_size, &type);
	//reconstruct_HuffTree_and_Decode_16states(tdps->typeArray, dataSeriesLength, &type);
	//memcpy(type, tdps->typeArray, dataSeriesLength*sizeof(unsigned short));
	//type = tdps->typeArray;
	decode_withTree(tdps->typeArray, dataSeriesLength, type);
	
	unsigned char preBytes[8];
	unsigned char curBytes[8];
	
	memset(preBytes, 0, 8);

	int curByteIndex = 0;
	int reqBytesLength, resiBitsLength, resiBits; 
	unsigned char leadingNum;	
	double medianValue, exactData, predValue;
	
	reqBytesLength = tdps->reqLength/8;
	resiBitsLength = tdps->reqLength%8;
	medianValue = tdps->medianValue;
	
	int type_;
	for (i = 0; i < dataSeriesLength; i++) {
		type_ = type[i];
		switch (type_) {
		case 0:
			// compute resiBits
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 8);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}
			
			exactData = bytesToDouble(curBytes);
			(*data)[i] = exactData + medianValue;
			memcpy(preBytes,curBytes,8);
			break;
		default:
			predValue = 2 * (*data)[i-1] - (*data)[i-2];
			(*data)[i] = predValue + (type_-intvRadius)*interval;
			break;
		}
		//printf("%.30G\n",(*data)[i]);
	}
	free(leadNum);
	free(type);
	return;
}

void decompressDataSeries_double_2D(double** data, int r1, int r2, TightDataPointStorageD* tdps) 
{
	updateQuantizationInfo(tdps->intervals);
	//printf("tdps->intervals=%d, intvRadius=%d\n", tdps->intervals, intvRadius);
	
	int i, j, k = 0, p = 0, l = 0; // k is to track the location of residual_bit
	// in resiMidBits, p is to track the
	// byte_index of resiMidBits, l is for
	// leadNum
	int dataSeriesLength = r1*r2;
	//	printf ("%d %d\n", r1, r2);

	unsigned char* leadNum;
	double realPrecision = tdps->realPrecision;

	convertByteArray2IntArray_fast_2b(tdps->exactDataNum, tdps->leadNumArray, tdps->leadNumArray_size, &leadNum);

	*data = (double*)malloc(sizeof(double)*dataSeriesLength);

	int* type = (int*)malloc(dataSeriesLength*sizeof(int));
	//convertByteArray2IntArray_fast_3b(dataSeriesLength, tdps->typeArray, tdps->typeArray_size, &type);
	//reconstruct_HuffTree_and_Decode_16states(tdps->typeArray, dataSeriesLength, &type);
	//memcpy(type, tdps->typeArray, dataSeriesLength*sizeof(unsigned short));
	//type = tdps->typeArray;
	decode_withTree(tdps->typeArray, dataSeriesLength, type);

	unsigned char preBytes[8];
	unsigned char curBytes[8];

	memset(preBytes, 0, 8);

	int curByteIndex = 0;
	int reqBytesLength, resiBitsLength, resiBits; 
	unsigned char leadingNum;	
	double medianValue, exactData, predValue;
	int type_;

	reqBytesLength = tdps->reqLength/8;
	resiBitsLength = tdps->reqLength%8;
	medianValue = tdps->medianValue;

	double pred1D, pred2D;
	int ii, jj;

	/* Process Row-0, data 0 */

	// compute resiBits
	resiBits = 0;
	if (resiBitsLength != 0) {
		int kMod8 = k % 8;
		int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
		if (rightMovSteps > 0) {
			int code = getRightMovingCode(kMod8, resiBitsLength);
			resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
		} else if (rightMovSteps < 0) {
			int code1 = getLeftMovingCode(kMod8);
			int code2 = getRightMovingCode(kMod8, resiBitsLength);
			int leftMovSteps = -rightMovSteps;
			rightMovSteps = 8 - leftMovSteps;
			resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
			p++;
			resiBits = resiBits
					| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
		} else // rightMovSteps == 0
		{
			int code = getRightMovingCode(kMod8, resiBitsLength);
			resiBits = (tdps->residualMidBits[p] & code);
			p++;
		}
		k += resiBitsLength;
	}

	// recover the exact data
	memset(curBytes, 0, 8);
	leadingNum = leadNum[l++];
	memcpy(curBytes, preBytes, leadingNum);
	for (j = leadingNum; j < reqBytesLength; j++)
		curBytes[j] = tdps->exactMidBytes[curByteIndex++];
	if (resiBitsLength != 0) {
		unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
		curBytes[reqBytesLength] = resiByte;
	}

	exactData = bytesToDouble(curBytes);
	(*data)[0] = exactData + medianValue;
	memcpy(preBytes,curBytes,8);

	/* Process Row-0, data 1 */
	pred1D = (*data)[0];

	type_ = type[1]; 
	if (type_ != 0)
	{
		(*data)[1] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
	}
	else
	{
		// compute resiBits
		resiBits = 0;
		if (resiBitsLength != 0) {
			int kMod8 = k % 8;
			int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
			if (rightMovSteps > 0) {
				int code = getRightMovingCode(kMod8, resiBitsLength);
				resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
			} else if (rightMovSteps < 0) {
				int code1 = getLeftMovingCode(kMod8);
				int code2 = getRightMovingCode(kMod8, resiBitsLength);
				int leftMovSteps = -rightMovSteps;
				rightMovSteps = 8 - leftMovSteps;
				resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
				p++;
				resiBits = resiBits
						| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
			} else // rightMovSteps == 0
			{
				int code = getRightMovingCode(kMod8, resiBitsLength);
				resiBits = (tdps->residualMidBits[p] & code);
				p++;
			}
			k += resiBitsLength;
		}

		// recover the exact data
		memset(curBytes, 0, 8);
		leadingNum = leadNum[l++];
		memcpy(curBytes, preBytes, leadingNum);
		for (j = leadingNum; j < reqBytesLength; j++)
			curBytes[j] = tdps->exactMidBytes[curByteIndex++];
		if (resiBitsLength != 0) {
			unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
			curBytes[reqBytesLength] = resiByte;
		}
		
		exactData = bytesToDouble(curBytes);
		(*data)[1] = exactData + medianValue;
		memcpy(preBytes,curBytes,8);
	}

	/* Process Row-0, data 2 --> data r2-1 */
	for (jj = 2; jj < r2; jj++)
	{
		pred1D = 2*(*data)[jj-1] - (*data)[jj-2];

		type_ = type[jj];
		if (type_ != 0)
		{
			(*data)[jj] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else
		{
			// compute resiBits
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 8);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToDouble(curBytes);
			(*data)[jj] = exactData + medianValue;
			memcpy(preBytes,curBytes,8);
		}
	}

	int index;
	/* Process Row-1 --> Row-r1-1 */
	for (ii = 1; ii < r1; ii++)
	{
		/* Process row-ii data 0 */
		index = ii*r2;
		pred1D = (*data)[index-r2];

		type_ = type[index];
		if (type_ != 0)
		{
			(*data)[index] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else
		{
			// compute resiBits
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 8);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToDouble(curBytes);
			(*data)[index] = exactData + medianValue;
			memcpy(preBytes,curBytes,8);
		}

		/* Process row-ii data 1 --> r2-1*/
		for (jj = 1; jj < r2; jj++)
		{
			index = ii*r2+jj;
			pred2D = (*data)[index-1] + (*data)[index-r2] - (*data)[index-r2-1];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else
			{
				// compute resiBits
				resiBits = 0;
				if (resiBitsLength != 0) {
					int kMod8 = k % 8;
					int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
					if (rightMovSteps > 0) {
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
					} else if (rightMovSteps < 0) {
						int code1 = getLeftMovingCode(kMod8);
						int code2 = getRightMovingCode(kMod8, resiBitsLength);
						int leftMovSteps = -rightMovSteps;
						rightMovSteps = 8 - leftMovSteps;
						resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
						p++;
						resiBits = resiBits
								| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
					} else // rightMovSteps == 0
					{
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code);
						p++;
					}
					k += resiBitsLength;
				}

				// recover the exact data
				memset(curBytes, 0, 8);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToDouble(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,8);
			}
		}
	}

	free(leadNum);
	free(type);
	return;
}

void decompressDataSeries_double_3D(double** data, int r1, int r2, int r3, TightDataPointStorageD* tdps) 
{
	updateQuantizationInfo(tdps->intervals);
	int i, j, k = 0, p = 0, l = 0; // k is to track the location of residual_bit
	// in resiMidBits, p is to track the
	// byte_index of resiMidBits, l is for
	// leadNum
	int dataSeriesLength = r1*r2*r3;
//	printf ("%d %d %d\n", r1, r2, r3);

	unsigned char* leadNum;
	double realPrecision = tdps->realPrecision;

	convertByteArray2IntArray_fast_2b(tdps->exactDataNum, tdps->leadNumArray, tdps->leadNumArray_size, &leadNum);

	*data = (double*)malloc(sizeof(double)*dataSeriesLength);

	int* type = (int*)malloc(dataSeriesLength*sizeof(int));
	//convertByteArray2IntArray_fast_3b(dataSeriesLength, tdps->typeArray, tdps->typeArray_size, &type);
	//reconstruct_HuffTree_and_Decode_16states(tdps->typeArray, dataSeriesLength, &type);
	//memcpy(type, tdps->typeArray, dataSeriesLength*sizeof(unsigned short));
	//type = tdps->typeArray;
	decode_withTree(tdps->typeArray, dataSeriesLength, type);
	//for(i=0;i<dataSeriesLength;i++)
	//	printf("%u\n", type[i]);

	unsigned char preBytes[8];
	unsigned char curBytes[8];

	memset(preBytes, 0, 8);

	int curByteIndex = 0;
	int reqBytesLength, resiBitsLength, resiBits;
	unsigned char leadingNum;
	double medianValue, exactData, predValue;
	int type_;

	reqBytesLength = tdps->reqLength/8;
	resiBitsLength = tdps->reqLength%8;
	medianValue = tdps->medianValue;

	double pred1D, pred2D, pred3D;
	int ii, jj, kk;

	///////////////////////////	Process layer-0 ///////////////////////////
	/* Process Row-0 data 0*/

	// compute resiBits
	resiBits = 0;
	if (resiBitsLength != 0) {
		int kMod8 = k % 8;
		int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
		if (rightMovSteps > 0) {
			int code = getRightMovingCode(kMod8, resiBitsLength);
			resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
		} else if (rightMovSteps < 0) {
			int code1 = getLeftMovingCode(kMod8);
			int code2 = getRightMovingCode(kMod8, resiBitsLength);
			int leftMovSteps = -rightMovSteps;
			rightMovSteps = 8 - leftMovSteps;
			resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
			p++;
			resiBits = resiBits
					| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
		} else // rightMovSteps == 0
		{
			int code = getRightMovingCode(kMod8, resiBitsLength);
			resiBits = (tdps->residualMidBits[p] & code);
			p++;
		}
		k += resiBitsLength;
	}

	// recover the exact data
	memset(curBytes, 0, 8);
	leadingNum = leadNum[l++];
	memcpy(curBytes, preBytes, leadingNum);
	for (j = leadingNum; j < reqBytesLength; j++)
		curBytes[j] = tdps->exactMidBytes[curByteIndex++];
	if (resiBitsLength != 0) {
		unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
		curBytes[reqBytesLength] = resiByte;
	}

	exactData = bytesToDouble(curBytes);
	(*data)[0] = exactData + medianValue;
	memcpy(preBytes,curBytes,8);

	/* Process Row-0, data 1 */
	pred1D = (*data)[0];

	type_ = type[1];
	if (type_ != 0)
	{
		(*data)[1] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
	}
	else
	{
		// compute resiBits
		resiBits = 0;
		if (resiBitsLength != 0) {
			int kMod8 = k % 8;
			int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
			if (rightMovSteps > 0) {
				int code = getRightMovingCode(kMod8, resiBitsLength);
				resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
			} else if (rightMovSteps < 0) {
				int code1 = getLeftMovingCode(kMod8);
				int code2 = getRightMovingCode(kMod8, resiBitsLength);
				int leftMovSteps = -rightMovSteps;
				rightMovSteps = 8 - leftMovSteps;
				resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
				p++;
				resiBits = resiBits
						| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
			} else // rightMovSteps == 0
			{
				int code = getRightMovingCode(kMod8, resiBitsLength);
				resiBits = (tdps->residualMidBits[p] & code);
				p++;
			}
			k += resiBitsLength;
		}

		// recover the exact data
		memset(curBytes, 0, 8);
		leadingNum = leadNum[l++];
		memcpy(curBytes, preBytes, leadingNum);
		for (j = leadingNum; j < reqBytesLength; j++)
			curBytes[j] = tdps->exactMidBytes[curByteIndex++];
		if (resiBitsLength != 0) {
			unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
			curBytes[reqBytesLength] = resiByte;
		}

		exactData = bytesToDouble(curBytes);
		(*data)[1] = exactData + medianValue;
		memcpy(preBytes,curBytes,8);
	}

	/* Process Row-0, data 2 --> data r3-1 */
	for (jj = 2; jj < r3; jj++)
	{
		pred1D = 2*(*data)[jj-1] - (*data)[jj-2];

		type_ = type[jj];
		if (type_ != 0)
		{
			(*data)[jj] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else
		{
			// compute resiBits
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 8);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToDouble(curBytes);
			(*data)[jj] = exactData + medianValue;
			memcpy(preBytes,curBytes,8);
		}
	}

	int index;
	/* Process Row-1 --> Row-r2-1 */
	for (ii = 1; ii < r2; ii++)
	{
		/* Process row-ii data 0 */
		index = ii*r3;
		pred1D = (*data)[index-r3];

		type_ = type[index];
		if (type_ != 0)
		{
			(*data)[index] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else
		{
			// compute resiBits
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 8);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToDouble(curBytes);
			(*data)[index] = exactData + medianValue;
			memcpy(preBytes,curBytes,8);
		}

		/* Process row-ii data 1 --> r3-1*/
		for (jj = 1; jj < r3; jj++)
		{
			index = ii*r3+jj;
			pred2D = (*data)[index-1] + (*data)[index-r3] - (*data)[index-r3-1];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else
			{
				// compute resiBits
				resiBits = 0;
				if (resiBitsLength != 0) {
					int kMod8 = k % 8;
					int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
					if (rightMovSteps > 0) {
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
					} else if (rightMovSteps < 0) {
						int code1 = getLeftMovingCode(kMod8);
						int code2 = getRightMovingCode(kMod8, resiBitsLength);
						int leftMovSteps = -rightMovSteps;
						rightMovSteps = 8 - leftMovSteps;
						resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
						p++;
						resiBits = resiBits
								| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
					} else // rightMovSteps == 0
					{
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code);
						p++;
					}
					k += resiBitsLength;
				}

				// recover the exact data
				memset(curBytes, 0, 8);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToDouble(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,8);
			}
		}
	}

	///////////////////////////	Process layer-1 --> layer-r1-1 ///////////////////////////

	for (kk = 1; kk < r1; kk++)
	{
		/* Process Row-0 data 0*/
		index = kk*r2*r3;
		pred1D = (*data)[index-r2*r3];

		type_ = type[index];
		if (type_ != 0)
		{
			(*data)[index] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else
		{
			// compute resiBits
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 8);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToDouble(curBytes);
			(*data)[index] = exactData + medianValue;
			memcpy(preBytes,curBytes,8);
		}

		/* Process Row-0 data 1 --> data r3-1 */
		for (jj = 1; jj < r3; jj++)
		{
			index = kk*r2*r3+jj;
			pred2D = (*data)[index-1] + (*data)[index-r2*r3] - (*data)[index-r2*r3-1];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else
			{
				// compute resiBits
				resiBits = 0;
				if (resiBitsLength != 0) {
					int kMod8 = k % 8;
					int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
					if (rightMovSteps > 0) {
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
					} else if (rightMovSteps < 0) {
						int code1 = getLeftMovingCode(kMod8);
						int code2 = getRightMovingCode(kMod8, resiBitsLength);
						int leftMovSteps = -rightMovSteps;
						rightMovSteps = 8 - leftMovSteps;
						resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
						p++;
						resiBits = resiBits
								| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
					} else // rightMovSteps == 0
					{
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code);
						p++;
					}
					k += resiBitsLength;
				}

				// recover the exact data
				memset(curBytes, 0, 8);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToDouble(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,8);
			}
		}

		/* Process Row-1 --> Row-r2-1 */
		for (ii = 1; ii < r2; ii++)
		{
			/* Process Row-i data 0 */
			index = kk*r2*r3 + ii*r3;
			pred2D = (*data)[index-r3] + (*data)[index-r2*r3] - (*data)[index-r2*r3-r3];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else
			{
				// compute resiBits
				resiBits = 0;
				if (resiBitsLength != 0) {
					int kMod8 = k % 8;
					int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
					if (rightMovSteps > 0) {
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
					} else if (rightMovSteps < 0) {
						int code1 = getLeftMovingCode(kMod8);
						int code2 = getRightMovingCode(kMod8, resiBitsLength);
						int leftMovSteps = -rightMovSteps;
						rightMovSteps = 8 - leftMovSteps;
						resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
						p++;
						resiBits = resiBits
								| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
					} else // rightMovSteps == 0
					{
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code);
						p++;
					}
					k += resiBitsLength;
				}

				// recover the exact data
				memset(curBytes, 0, 8);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToDouble(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,8);
			}

			/* Process Row-i data 1 --> data r3-1 */
			for (jj = 1; jj < r3; jj++)
			{
				index = kk*r2*r3 + ii*r3 + jj;
				pred3D = (*data)[index-1] + (*data)[index-r3] + (*data)[index-r2*r3]
					- (*data)[index-r3-1] - (*data)[index-r2*r3-r3] - (*data)[index-r2*r3-1] + (*data)[index-r2*r3-r3-1];

				type_ = type[index];
				if (type_ != 0)
				{
					(*data)[index] = pred3D + 2 * (type_ - intvRadius) * realPrecision;
				}
				else
				{
					// compute resiBits
					resiBits = 0;
					if (resiBitsLength != 0) {
						int kMod8 = k % 8;
						int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
						if (rightMovSteps > 0) {
							int code = getRightMovingCode(kMod8, resiBitsLength);
							resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
						} else if (rightMovSteps < 0) {
							int code1 = getLeftMovingCode(kMod8);
							int code2 = getRightMovingCode(kMod8, resiBitsLength);
							int leftMovSteps = -rightMovSteps;
							rightMovSteps = 8 - leftMovSteps;
							resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
							p++;
							resiBits = resiBits
									| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
						} else // rightMovSteps == 0
						{
							int code = getRightMovingCode(kMod8, resiBitsLength);
							resiBits = (tdps->residualMidBits[p] & code);
							p++;
						}
						k += resiBitsLength;
					}

					// recover the exact data
					memset(curBytes, 0, 8);
					leadingNum = leadNum[l++];
					memcpy(curBytes, preBytes, leadingNum);
					for (j = leadingNum; j < reqBytesLength; j++)
						curBytes[j] = tdps->exactMidBytes[curByteIndex++];
					if (resiBitsLength != 0) {
						unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
						curBytes[reqBytesLength] = resiByte;
					}

					exactData = bytesToDouble(curBytes);
					(*data)[index] = exactData + medianValue;
					memcpy(preBytes,curBytes,8);
				}
			}
		}
	}

	free(leadNum);
	free(type);
	return;
}

void getSnapshotData_double_1D(double** data, int dataSeriesLength, TightDataPointStorageD* tdps) 
{	
	SZ_Reset();
	int i;
	if (tdps->allSameData) {
		double value = bytesToDouble(tdps->exactMidBytes);
		*data = (double*)malloc(sizeof(double)*dataSeriesLength);
		for (i = 0; i < dataSeriesLength; i++)
			(*data)[i] = value;
	} else {
		if (tdps->rtypeArray == NULL) {
			
			decompressDataSeries_double_1D(data, dataSeriesLength, tdps);
			return;
		} else {
			*data = (double*)malloc(sizeof(double)*dataSeriesLength);
			// insert the reserved values
			//int[] rtypes = TypeManager.convertByteArray2IntArray_fast_1b(
			//		dataSeriesLength, rtypeArray);
			int* rtypes;
			int validLength = computeBitNumRequired(dataSeriesLength);
			decompressBitArraybySimpleLZ77(&rtypes, tdps->rtypeArray, tdps->rtypeArray_size, dataSeriesLength, validLength);
			int count = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 1)
					(*data)[i] = tdps->reservedValue;
				else
					count++;
			}
			// get the decompressed data
			double* decmpData;
			decompressDataSeries_double_1D(&decmpData, dataSeriesLength, tdps);
			// insert the decompressed data
			int k = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 0) {
					(*data)[i] = decmpData[k++];
				}
			}
			free(decmpData);
			free(rtypes);
		}
	}
}

void getSnapshotData_double_2D(double** data, int r1, int r2, TightDataPointStorageD* tdps) 
{
	SZ_Reset();
	int i;
	int dataSeriesLength = r1*r2;
	if (tdps->allSameData) {
		double value = bytesToDouble(tdps->exactMidBytes);
		*data = (double*)malloc(sizeof(double)*dataSeriesLength);
		for (i = 0; i < dataSeriesLength; i++)
			(*data)[i] = value;
	} else {
		if (tdps->rtypeArray == NULL) {

			decompressDataSeries_double_2D(data, r1, r2, tdps);
			return;
		} else {
			*data = (double*)malloc(sizeof(double)*dataSeriesLength);
			// insert the reserved values
			//int[] rtypes = TypeManager.convertByteArray2IntArray_fast_1b(
			//		dataSeriesLength, rtypeArray);
			int* rtypes;
			int validLength = computeBitNumRequired(dataSeriesLength);
			decompressBitArraybySimpleLZ77(&rtypes, tdps->rtypeArray, tdps->rtypeArray_size, dataSeriesLength, validLength);
			int count = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 1)
					(*data)[i] = tdps->reservedValue;
				else
					count++;
			}
			// get the decompressed data
			double* decmpData;
			decompressDataSeries_double_2D(&decmpData, r1, r2, tdps);
			// insert the decompressed data
			int k = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 0) {
					(*data)[i] = decmpData[k++];
				}
			}
			free(decmpData);
			free(rtypes);
		}
	}
}

void getSnapshotData_double_3D(double** data, int r1, int r2, int r3, TightDataPointStorageD* tdps) 
{
	SZ_Reset();
	int i;
	int dataSeriesLength = r1*r2*r3;
	if (tdps->allSameData) {
		double value = bytesToDouble(tdps->exactMidBytes);
		*data = (double*)malloc(sizeof(double)*dataSeriesLength);
		for (i = 0; i < dataSeriesLength; i++)
			(*data)[i] = value;
	} else {
		if (tdps->rtypeArray == NULL) {

			decompressDataSeries_double_3D(data, r1, r2, r3, tdps);
			return;
		} else {
			*data = (double*)malloc(sizeof(double)*dataSeriesLength);
			// insert the reserved values
			//int[] rtypes = TypeManager.convertByteArray2IntArray_fast_1b(
			//		dataSeriesLength, rtypeArray);
			int* rtypes;
			int validLength = computeBitNumRequired(dataSeriesLength);
			decompressBitArraybySimpleLZ77(&rtypes, tdps->rtypeArray, tdps->rtypeArray_size, dataSeriesLength, validLength);
			int count = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 1)
					(*data)[i] = tdps->reservedValue;
				else
					count++;
			}
			// get the decompressed data
			double* decmpData;
			decompressDataSeries_double_3D(&decmpData, r1, r2, r3, tdps);
			// insert the decompressed data
			int k = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 0) {
					(*data)[i] = decmpData[k++];
				}
			}
			free(decmpData);
			free(rtypes);
		}
	}
}

/**
 * 
 * type's length == dataSeriesLength
 * exactMidBytes's length == exactMidBytes_size
 * leadNumIntArray's length == exactDataNum
 * escBytes's length == escBytes_size
 * resiBitLength's length == resiBitLengthSize
 * */
void new_TightDataPointStorageD(TightDataPointStorageD **this, 
		int dataSeriesLength, int exactDataNum, 
		int* type, unsigned char* exactMidBytes, int exactMidBytes_size,
		unsigned char* leadNumIntArray,  //leadNumIntArray contains readable numbers....
		unsigned char* resiMidBits, int resiMidBits_size,
		unsigned char* resiBitLength, int resiBitLengthSize, 
		double realPrecision, double medianValue, char reqLength, unsigned int intervals) {
	int i = 0;
	*this = (TightDataPointStorageD *)malloc(sizeof(TightDataPointStorageD));
	(*this)->allSameData = 0;
	(*this)->realPrecision = realPrecision;
	(*this)->medianValue = medianValue;
	(*this)->reqLength = reqLength;
	
	(*this)->dataSeriesLength = dataSeriesLength;
	(*this)->exactDataNum = exactDataNum;
	
	(*this)->rtypeArray = NULL;
	(*this)->rtypeArray_size = 0;
	
	//(*this)->typeArray_size = convertIntArray2ByteArray_fast_3b(type, dataSeriesLength, &((*this)->typeArray));
	//for(;i<dataSeriesLength;i++)
	//	type[i]+=48;
	//huff_init(type);
	
	//(*this)->typeArray_size = convert_HuffTree_and_Encode_16states(type, dataSeriesLength, &((*this)->typeArray));
	
	/*(*this)->typeArray_size = sizeof(unsigned short)*dataSeriesLength;
	(*this)->typeArray = (unsigned char*)malloc(dataSeriesLength*sizeof(unsigned short));
	memcpy((*this)->typeArray, type, dataSeriesLength*sizeof(unsigned short));
	free(type);*/
	
	//for(i=0;i<dataSeriesLength;i++)
	//	printf("%u\n", type[i]);
	encode_withTree(type, dataSeriesLength, &(*this)->typeArray, &(*this)->typeArray_size);

	(*this)->exactMidBytes = exactMidBytes;
	(*this)->exactMidBytes_size = exactMidBytes_size;
	
	(*this)->leadNumArray_size = convertIntArray2ByteArray_fast_2b(leadNumIntArray, exactDataNum, &((*this)->leadNumArray));

	//(*this)->residualMidBits = resiMidBits;
	//(*this)->residualMidBits_size = resiMidBits_size;

	(*this)->residualMidBits_size = convertIntArray2ByteArray_fast_dynamic(resiMidBits, resiBitLength, resiBitLengthSize, &((*this)->residualMidBits));
	
	(*this)->intervals = intervals;
	(*this)->isLossless = 0;
}

//TODO: convert TightDataPointStorageD to bytes...
void convertTDPStoFlatBytes_double(TightDataPointStorageD *tdps, unsigned char** bytes, int *size) 
{
	int i, k = 0; 
	unsigned char intervalsBytes[4];
	unsigned char typeArrayLengthBytes[4];
	unsigned char rTypeLengthBytes[4];
	unsigned char dsLengthBytes[4];
	unsigned char exactLengthBytes[4];
	unsigned char exactMidBytesLength[4];
	unsigned char reservedValueBytes[8];
	unsigned char realPrecisionBytes[8];
	
	unsigned char medianValueBytes[8];

	intToBytes_bigEndian(dsLengthBytes, tdps->dataSeriesLength);//4
	unsigned char sameByte = tdps->allSameData==1?(unsigned char)1:(unsigned char)0;
	sameByte = sameByte | (szMode << 1);
	
	if(tdps->allSameData==1)
	{
		int totalByteLength = 3 + 4 + 1 + tdps->exactMidBytes_size;
		*bytes = (unsigned char *)malloc(sizeof(unsigned char)*totalByteLength);
	
		for (i = 0; i < 3; i++)//3
			(*bytes)[k++] = versionNumber[i];
		for (i = 0; i < 4; i++)
			(*bytes)[k++] = dsLengthBytes[i];
		(*bytes)[k++] = sameByte;
		for (i = 0; i < tdps->exactMidBytes_size; i++)
			(*bytes)[k++] = tdps->exactMidBytes[i];
		
		*size = totalByteLength;
	}
	else if (tdps->rtypeArray == NULL) 
	{
		int residualMidBitsLength = tdps->residualMidBits == NULL ? 0 : tdps->residualMidBits_size;
		int totalByteLength = 3 + 4 + 1 + 4 + 8 + 1 + 8 + 4 + 4 + 4 
				+ tdps->typeArray_size + tdps->leadNumArray_size + tdps->exactMidBytes_size + residualMidBitsLength;

		*bytes = (unsigned char *)malloc(sizeof(unsigned char)*totalByteLength);

		for(i = 0;i<3;i++)//3 bytes
			(*bytes)[k++] = versionNumber[i];
		for(i = 0;i<4;i++)//4 bytes
			(*bytes)[k++] = dsLengthBytes[i];
			
		(*bytes)[k++] = sameByte;	//1	byte		
		
		intToBytes_bigEndian(intervalsBytes, tdps->intervals);
		for(i = 0;i<4;i++)//4
			(*bytes)[k++] = intervalsBytes[i];		
		
		doubleToBytes(medianValueBytes, tdps->medianValue);
		for (i = 0; i < 8; i++)// 8
			(*bytes)[k++] = medianValueBytes[i];		

		(*bytes)[k++] = tdps->reqLength; //1 byte

		doubleToBytes(realPrecisionBytes, tdps->realPrecision);
		for (i = 0; i < 8; i++)// 8
			(*bytes)[k++] = realPrecisionBytes[i];
				
		intToBytes_bigEndian(typeArrayLengthBytes, tdps->typeArray_size);
		for(i = 0;i<4;i++)//4
			(*bytes)[k++] = typeArrayLengthBytes[i];				
					
		intToBytes_bigEndian(exactLengthBytes, tdps->exactDataNum);
		for(i = 0;i<4;i++)//4
			(*bytes)[k++] = exactLengthBytes[i];

		intToBytes_bigEndian(exactMidBytesLength, tdps->exactMidBytes_size);
		for(i = 0;i<4;i++)//4
			(*bytes)[k++] = exactMidBytesLength[i];

		memcpy(&((*bytes)[k]), tdps->typeArray, tdps->typeArray_size);
		k += tdps->typeArray_size;
		memcpy(&((*bytes)[k]), tdps->leadNumArray, tdps->leadNumArray_size);
		k += tdps->leadNumArray_size;
		memcpy(&((*bytes)[k]), tdps->exactMidBytes, tdps->exactMidBytes_size);
		k += tdps->exactMidBytes_size;		
		if(tdps->residualMidBits!=NULL)
		{
			memcpy(&((*bytes)[k]), tdps->residualMidBits, tdps->residualMidBits_size);
			k += tdps->residualMidBits_size;	
		}

		*size = totalByteLength;
	}
	else //the case with reserved value
	{
		int residualMidBitsLength = tdps->residualMidBits == NULL ? 0 : tdps->residualMidBits_size;
		int totalByteLength = 3 + 4 + 1 + 4 + 8 + 1 + 8 + 4 + 4 + 4 + 4 + 8 + tdps->rtypeArray_size
		+ tdps->typeArray_size + tdps->leadNumArray_size 
		+ tdps->exactMidBytes_size + residualMidBitsLength;

		sameByte = (unsigned char) (sameByte | 0x08); // 00001000, the 4th bit
												// denotes whether it is
												// with "reserved value"

		*bytes = (unsigned char*)malloc(sizeof(unsigned char)*totalByteLength);
		
		for(i = 0;i<3;i++)//3
			(*bytes)[k++] = versionNumber[i];
		for(i = 0;i<4;i++)//4
			(*bytes)[k++] = dsLengthBytes[i];		
			
		(*bytes)[k++] = sameByte;						//1

		intToBytes_bigEndian(intervalsBytes, tdps->intervals);
		for(i = 0;i<4;i++)//4
			(*bytes)[k++] = intervalsBytes[i];	

		doubleToBytes(medianValueBytes, tdps->medianValue);
		for (i = 0; i < 8; i++)// 8
			(*bytes)[k++] = medianValueBytes[i];		

		(*bytes)[k++] = tdps->reqLength; //1 byte

		doubleToBytes(realPrecisionBytes, tdps->realPrecision);
		for (i = 0; i < 8; i++)// 8
			(*bytes)[k++] = realPrecisionBytes[i];		
		
		intToBytes_bigEndian(typeArrayLengthBytes, tdps->typeArray_size);
		for(i = 0;i<4;i++)//4
			(*bytes)[k++] = typeArrayLengthBytes[i];			
		
		intToBytes_bigEndian(rTypeLengthBytes, tdps->rtypeArray_size);
		for(i = 0;i<4;i++)//4
			(*bytes)[k++] = rTypeLengthBytes[i];	
		
		intToBytes_bigEndian(exactLengthBytes, tdps->exactDataNum);
		for(i = 0;i<4;i++)//4
			(*bytes)[k++] = exactLengthBytes[i];

		intToBytes_bigEndian(exactMidBytesLength, tdps->exactMidBytes_size);
		for(i = 0;i<4;i++)//4
			(*bytes)[k++] = exactMidBytesLength[i];

		doubleToBytes(reservedValueBytes, tdps->reservedValue);
		for (i = 0; i < 8; i++)// 8
			(*bytes)[k++] = reservedValueBytes[i];
		
		memcpy(&((*bytes)[k]), tdps->rtypeArray, tdps->rtypeArray_size);
		k += tdps->rtypeArray_size;		
		memcpy(&((*bytes)[k]), tdps->typeArray, tdps->typeArray_size);
		k += tdps->typeArray_size;
		memcpy(&((*bytes)[k]), tdps->leadNumArray, tdps->leadNumArray_size);
		k += tdps->leadNumArray_size;
		memcpy(&((*bytes)[k]), tdps->exactMidBytes, tdps->exactMidBytes_size);
		k += tdps->exactMidBytes_size;		
		if(tdps->residualMidBits!=NULL)
		{
			memcpy(&((*bytes)[k]), tdps->residualMidBits, tdps->residualMidBits_size);
			k += tdps->residualMidBits_size;	
		}

		*size = totalByteLength;
	}
}

void free_TightDataPointStorageD(TightDataPointStorageD *tdps)
{
	if(tdps->rtypeArray!=NULL)
		free(tdps->rtypeArray);
	if(tdps->typeArray!=NULL)
		free(tdps->typeArray);
	if(tdps->leadNumArray!=NULL)
		free(tdps->leadNumArray);
//	if(tdps->exactMidBytes!=NULL)
//		free(tdps->exactMidBytes);
	//if(tdps->escBytes!=NULL)
	//	free(tdps->escBytes);
	if(tdps->residualMidBits!=NULL)
		free(tdps->residualMidBits);
	free(tdps);
}
