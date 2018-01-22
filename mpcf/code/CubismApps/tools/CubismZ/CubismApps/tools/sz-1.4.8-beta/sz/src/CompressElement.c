/**
 *  @file CompressElement.c
 *  @author Sheng Di
 *  @date May, 2016
 *  @brief Functions of CompressElement
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdlib.h> 
#include <stdio.h>
#include <sz.h>
#include <CompressElement.h>

inline void listAdd_double(double last3CmprsData[3], double value)
{
	last3CmprsData[2] = last3CmprsData[1];
	last3CmprsData[1] = last3CmprsData[0];
	last3CmprsData[0] = value;
}

inline void listAdd_float(float last3CmprsData[3], float value)
{
	last3CmprsData[2] = last3CmprsData[1];
	last3CmprsData[1] = last3CmprsData[0];
	last3CmprsData[0] = value;
}

int validPrediction_double(double minErr, double precision)
{
	if(minErr<=precision)
		return 1;
	else
		return 0;
}

int validPrediction_float(float minErr, float precision)
{
	if(minErr<=precision)
		return 1;
	else
		return 0;
}

void new_LossyCompressionElement(LossyCompressionElement *lce, int leadingNum, unsigned char* intMidBytes, 
int intMidBytes_Length, int resiMidBitsLength, int resiBits)
{
	lce->leadingZeroBytes = leadingNum; //0,1,2,or 3
	memcpy(lce->integerMidBytes,intMidBytes,intMidBytes_Length);
	lce->integerMidBytes_Length = intMidBytes_Length; //they are mid_bits actually
	lce->resMidBitsLength = resiMidBitsLength;
	lce->residualMidBits = resiBits;
}

void updateLossyCompElement_Double(unsigned char* curBytes, unsigned char* preBytes, 
		int reqBytesLength, int resiBitsLength,  LossyCompressionElement *lce)
{
	int i, resiIndex, intMidBytes_Length = 0;
	int leadingNum = compIdenticalLeadingBytesCount_double(preBytes, curBytes); //in fact, float is enough for both single-precision and double-precisiond ata.
	int fromByteIndex = leadingNum;
	int toByteIndex = reqBytesLength; //later on: should use "< toByteIndex" to tarverse....
	if(fromByteIndex < toByteIndex)
	{
		intMidBytes_Length = reqBytesLength - leadingNum;
		memcpy(lce->integerMidBytes, &(curBytes[fromByteIndex]), intMidBytes_Length);
	}
	int resiBits = 0;
	if(resiBitsLength!=0)
	{
		resiIndex = reqBytesLength;
		if(resiIndex < 8)
			resiBits = (curBytes[resiIndex] & 0xFF) >> (8-resiBitsLength);
	}
	lce->leadingZeroBytes = leadingNum;
	lce->integerMidBytes_Length = intMidBytes_Length;
	lce->resMidBitsLength = resiBitsLength;
	lce->residualMidBits = resiBits;
}

void updateLossyCompElement_Float(unsigned char* curBytes, unsigned char* preBytes, 
		int reqBytesLength, int resiBitsLength,  LossyCompressionElement *lce)
{
	int i, resiIndex, intMidBytes_Length = 0;
	int leadingNum = compIdenticalLeadingBytesCount_float(preBytes, curBytes); //in fact, float is enough for both single-precision and double-precisiond ata.
	int fromByteIndex = leadingNum;
	int toByteIndex = reqBytesLength; //later on: should use "< toByteIndex" to tarverse....
	if(fromByteIndex < toByteIndex)
	{
		intMidBytes_Length = reqBytesLength - leadingNum;
		memcpy(lce->integerMidBytes, &(curBytes[fromByteIndex]), intMidBytes_Length);
	}
	int resiBits = 0;
	if(resiBitsLength!=0)
	{
		resiIndex = reqBytesLength;
		if(resiIndex < 8)
			resiBits = (curBytes[resiIndex] & 0xFF) >> (8-resiBitsLength);
	}
	lce->leadingZeroBytes = leadingNum;
	lce->integerMidBytes_Length = intMidBytes_Length;
	lce->resMidBitsLength = resiBitsLength;
	lce->residualMidBits = resiBits;
}
