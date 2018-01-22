/**
 *  @file TypeManager.c
 *  @author Sheng Di
 *  @date May, 2016
 *  @brief TypeManager is used to manage the type array: parsing of the bytes and other types in between.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include "DynamicByteArray.h"
#include "sz.h"

//int convertIntArray2ByteArray_fast_8b()

/**
 * little endian
 * [01|10|11|00|....]-->[01|10|11|00][....]
 * @param timeStepType
 * @return
 */
int convertIntArray2ByteArray_fast_2b(unsigned char* timeStepType, int timeStepTypeLength, unsigned char **result)
{
	int i, j, byteLength = 0;
	if(timeStepTypeLength%4==0)
		byteLength = timeStepTypeLength*2/8;
	else
		byteLength = timeStepTypeLength*2/8+1;
	if(byteLength>0)
		*result = (unsigned char*)malloc(byteLength*sizeof(unsigned char));
	else
		*result = NULL;
	int n = 0;
	for(i = 0;i<byteLength;i++)
	{
		int tmp = 0;
		for(j = 0;j<4&&n<timeStepTypeLength;j++)
		{
			int type = timeStepType[n];
			switch(type)
			{
			case 0: 
				
				break;
			case 1:
				tmp = (tmp | (1 << (6-j*2)));
				break;
			case 2:
				tmp = (tmp | (2 << (6-j*2)));
				break;
			case 3:
				tmp = (tmp | (3 << (6-j*2)));
				break;
			default:
				printf("Error: wrong timestep type...: type[%d]=%d\n", n, type);
				exit(0);
			}
			n++;
		}
		(*result)[i] = (unsigned char)tmp;
	}
	return byteLength;
}

void convertByteArray2IntArray_fast_2b(int stepLength, unsigned char* byteArray, int byteArrayLength, unsigned char **intArray)
{
	if(stepLength > byteArrayLength*4)
	{
		printf("Error: stepLength > byteArray.length*4\n");
		printf("stepLength=%d, byteArray.length=%d\n", stepLength, byteArrayLength);
		exit(0);
	}
	if(stepLength>0)
		*intArray = (unsigned char*)malloc(stepLength*sizeof(unsigned char));
	else
		*intArray = NULL;
	int i, n = 0;

	for (i = 0; i < byteArrayLength; i++) {
		unsigned char tmp = byteArray[i];
		(*intArray)[n++] = (tmp & 0xC0) >> 6;
		if(n==stepLength)
			break;
		(*intArray)[n++] = (tmp & 0x30) >> 4;
		if(n==stepLength)
			break;
		(*intArray)[n++] = (tmp & 0x0C) >> 2;
		if(n==stepLength)
			break;
		(*intArray)[n++] = tmp & 0x03;
		if(n==stepLength)
			break;
	}
}

int convertIntArray2ByteArray_fast_3b(unsigned char* timeStepType, int timeStepTypeLength, unsigned char **result)
{	
	int i = 0, j = 0, k = 0, byteLength = 0;
	if(timeStepTypeLength%8==0)
		byteLength = timeStepTypeLength*3/8;
	else
		byteLength = timeStepTypeLength*3/8+1;

	if(byteLength>0)
		*result = (unsigned char*)malloc(byteLength*sizeof(unsigned char));
	else
		*result = NULL;
	int n = 0, tmp = 0;
	for(n = 0;n<timeStepTypeLength;n++)
	{
		k = n%8;
		switch(k)
		{
		case 0:
			tmp = tmp | (timeStepType[n] << 5);
			break;
		case 1:
			tmp = tmp | (timeStepType[n] << 2);
			break;
		case 2: 
			tmp = tmp | (timeStepType[n] >> 1);
			(*result)[i++] = (unsigned char)tmp;
			tmp = 0 | (timeStepType[n] << 7);
			break;
		case 3:
			tmp = tmp | (timeStepType[n] << 4);
			break;
		case 4:
			tmp = tmp | (timeStepType[n] << 1);
			break;
		case 5:
			tmp = tmp | (timeStepType[n] >> 2);
			(*result)[i++] = (unsigned char)tmp;
			tmp = 0 | (timeStepType[n] << 6);
			break;
		case 6:
			tmp = tmp | (timeStepType[n] << 3);
			break;
		case 7:
			tmp = tmp | (timeStepType[n] << 0);
			(*result)[i++] = (unsigned char)tmp;
			tmp = 0;
			break;
		}
	}
	if(k!=7) //load the last one
		(*result)[i] = (unsigned char)tmp;
	
	return byteLength;
}

void convertByteArray2IntArray_fast_3b(int stepLength, unsigned char* byteArray, int byteArrayLength, unsigned char **intArray)
{	
	if(stepLength > byteArrayLength*8/3)
	{
		printf("Error: stepLength > byteArray.length*8/3\n");
		printf("stepLength=%d, byteArray.length=%d\n", stepLength, byteArrayLength);
		exit(0);		
	}
	if(stepLength>0)
		*intArray = (unsigned char*)malloc(stepLength*sizeof(unsigned char));
	else
		*intArray = NULL;
	int i = 0, ii = 0, n = 0;
	unsigned char tmp = byteArray[i];	
	for(n=0;n<stepLength;)
	{
		switch(n%8)
		{
		case 0:
			(*intArray)[n++] = (tmp & 0xE0) >> 5;
			break;
		case 1: 
			(*intArray)[n++] = (tmp & 0x1C) >> 2;
			break;
		case 2:
			ii = (tmp & 0x03) << 1;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0x80) >> 7;
			(*intArray)[n++] = ii;
			break;
		case 3:
			(*intArray)[n++] = (tmp & 0x70) >> 4;
			break;
		case 4:
			(*intArray)[n++] = (tmp & 0x0E) >> 1;
			break;
		case 5:
			ii = (tmp & 0x01) << 2;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0xC0) >> 6;
			(*intArray)[n++] = ii;
			break;
		case 6: 
			(*intArray)[n++] = (tmp & 0x38) >> 3;
			break;
		case 7:
			(*intArray)[n++] = (tmp & 0x07);
			i++;
			tmp = byteArray[i];
			break;
		}
	}
}

int getLeftMovingSteps(int k, unsigned char resiBitLength)
{
	return 8 - k%8 - resiBitLength;
}

/**
 * 
 * @param timeStepType is the resiMidBits
 * @param resiBitLength is the length of resiMidBits for each element, (the number of resiBitLength == the # of unpredictable elements
 * @return
 */
int convertIntArray2ByteArray_fast_dynamic(unsigned char* timeStepType, unsigned char* resiBitLength, int resiBitLengthLength, unsigned char **bytes)
{
	int i = 0, j = 0, k = 0, value;
	DynamicByteArray* dba;
	new_DBA(&dba, 1024);
	int tmp = 0, leftMovSteps = 0;
	for(j = 0;j<resiBitLengthLength;j++)
	{
		unsigned char rbl = resiBitLength[j];
		if(rbl==0)
			continue;
		value = timeStepType[i];
		leftMovSteps = getLeftMovingSteps(k, rbl);
		if(leftMovSteps < 0)
		{
			tmp = tmp | (value >> (-leftMovSteps));
			addDBA_Data(dba, (unsigned char)tmp);
			tmp = 0 | (value << (8+leftMovSteps));
		}
		else if(leftMovSteps > 0)
		{
			tmp = tmp | (value << leftMovSteps);
		}
		else //==0
		{
			tmp = tmp | value;
			addDBA_Data(dba, (unsigned char)tmp);
			tmp = 0;
		}
		i++;
		k += rbl;
	}
	if(leftMovSteps != 0)
		addDBA_Data(dba, (unsigned char)tmp);
	convertDBAtoBytes(dba, bytes);
	int size = dba->size;
	free_DBA(dba);
	return size;
}

int computeBitNumRequired(int dataLength)
{
	return 32 - numberOfLeadingZeros_Int(dataLength);
}

void decompressBitArraybySimpleLZ77(int** result, unsigned char* bytes, int bytesLength, int totalLength, int validLength)
{
	int pairLength = (bytesLength*8)/(validLength+1);
	int tmpLength = pairLength*2;
	int tmpResult[tmpLength];
	int i, j, k = 0;
	for(i = 0;i<tmpLength;i+=2)
	{
		int outIndex = k/8;
		int innerIndex = k%8;

		unsigned char curByte = bytes[outIndex];
		tmpResult[i] = (curByte >> (8-1-innerIndex)) & 0x01;
		k++;
		
		int numResult = extractBytes(bytes, k, validLength);
		
		tmpResult[i+1] = numResult;
		k = k + validLength;
	}
	
	*result = (int*)malloc(sizeof(int)*totalLength);
	k = 0;
	for(i = 0;i<tmpLength;i=i+2)
	{
		int state = tmpResult[i];
		int num = tmpResult[i+1];
		for(j = 0;j<num;j++)
			(*result)[k++] = state;
	}
}
