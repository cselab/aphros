/**
 *  @file rw.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief io interface for fortrance
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rw.h"
#include "sz.h"

int checkFileSize(char *srcFilePath)
{
	int filesize;
	FILE *pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
	{
		printf("Failed to open input file. 1\n");
		exit(1);
	}
	fseek(pFile, 0, SEEK_END);
    filesize = (int)ftell(pFile);
    fclose(pFile);
    return filesize;
}

unsigned char *readByteData(char *srcFilePath, int *byteLength)
{
	FILE *pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 1\n");
        exit(1);
    }
	fseek(pFile, 0, SEEK_END);
    *byteLength = (int)ftell(pFile);
    fclose(pFile);
    
    unsigned char *byteBuf = ( unsigned char *)malloc((*byteLength)*sizeof(unsigned char)); //sizeof(char)==1
    
    pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 2\n");
        exit(1);
    }
    fread(byteBuf, 1, *byteLength, pFile);
    fclose(pFile);
    return byteBuf;
}

double *readDoubleData(char *srcFilePath, int *nbEle)
{
	if(dataEndianType==sysEndianType)
	{
		double *daBuf = readDoubleData_systemEndian(srcFilePath, nbEle);
		return daBuf;
	}
	else
	{
		int i,j;
		
		int byteLength;
		unsigned char* bytes = readByteData(srcFilePath, &byteLength);
		double *daBuf = (double *)malloc(byteLength);
		*nbEle = byteLength/8;
		
		ldouble buf;
		for(i = 0;i<*nbEle;i++)
		{
			j = i*8;
			memcpy(buf.byte, bytes+j, 8);
			symTransform_8bytes(buf.byte);
			daBuf[i] = buf.value;
		}
		free(bytes);
		return daBuf;
	}
}

float *readFloatData(char *srcFilePath, int *nbEle)
{
	if(dataEndianType==sysEndianType)
	{
		float *daBuf = readFloatData_systemEndian(srcFilePath, nbEle);
		return daBuf;
	}
	else
	{
		int i,j;
		
		int byteLength;
		unsigned char* bytes = readByteData(srcFilePath, &byteLength);
		float *daBuf = (float *)malloc(byteLength);
		*nbEle = byteLength/4;
		
		lfloat buf;
		for(i = 0;i<*nbEle;i++)
		{
			j = i*4;
			memcpy(buf.byte, bytes+j, 4);
			symTransform_4bytes(buf.byte);
			daBuf[i] = buf.value;
		}
		free(bytes);
		return daBuf;
	}
}

double *readDoubleData_systemEndian(char *srcFilePath, int *nbEle)
{
	int inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 1\n");
        exit(1);
    }
	fseek(pFile, 0, SEEK_END);
    inSize = ftell(pFile);
    *nbEle = (int)inSize/8; //only support double in this version
    fclose(pFile);
    
    double *daBuf = (double *)malloc(inSize);
    
    pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 2\n");
        exit(1);
    }
    fread(daBuf, 8, *nbEle, pFile);
    fclose(pFile);
    return daBuf;
}

float *readFloatData_systemEndian(char *srcFilePath, int *nbEle)
{
	int inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 1\n");
        exit(1);
    }
	fseek(pFile, 0, SEEK_END);
    inSize = ftell(pFile);
    *nbEle = (int)inSize/4; 
    fclose(pFile);
    
    if(inSize<=0)
    {
		printf("Error: input file is wrong!\n");
		exit(0);
	}
    
    float *daBuf = (float *)malloc(inSize);
    
    pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 2\n");
        exit(1);
    }
    fread(daBuf, 4, *nbEle, pFile);
    fclose(pFile);
    return daBuf;
}

void writeByteData(unsigned char *bytes, int byteLength, char *tgtFilePath)
{
	FILE *pFile = fopen(tgtFilePath, "wb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 3\n");
        exit(1);
    }
    
    fwrite(bytes, 1, byteLength, pFile); //write outSize bytes
    fclose(pFile);
}

void writeDoubleData(double *data, int nbEle, char *tgtFilePath)
{
	int i = 0;
	char s[64];
	FILE *pFile = fopen(tgtFilePath, "wb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 3\n");
        exit(1);
    }
    
    for(i = 0;i<nbEle;i++)
	{
		sprintf(s,"%.20G\n",data[i]);
		fputs(s, pFile);
	}
    
    fclose(pFile);
}

void writeFloatData(float *data, int nbEle, char *tgtFilePath)
{
	int i = 0;
	char s[64];
	FILE *pFile = fopen(tgtFilePath, "wb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 3\n");
        exit(1);
    }
   
    for(i = 0;i<nbEle;i++)
	{
		//printf("i=%d\n",i);
		//printf("data[i]=%f\n",data[i]);
		sprintf(s,"%.30G\n",data[i]);
		fputs(s, pFile);
	}
    
    fclose(pFile);
}

void writeData(void *data, int dataType, int nbEle, char *tgtFilePath)
{
	if(dataType == SZ_FLOAT)
	{
		float* dataArray = (float *)data;
		writeFloatData(dataArray, nbEle, tgtFilePath);
	}
	else if(dataType == SZ_DOUBLE)
	{
		double* dataArray = (double *)data;
		writeDoubleData(dataArray, nbEle, tgtFilePath);
	}
	else
	{
		printf("Error: data type cannot be the types other than SZ_FLOAT or SZ_DOUBLE\n");
		exit(0);	
	}
}

void writeFloatData_inBytes(float *data, int nbEle, char* tgtFilePath)
{
	int i = 0;
	lfloat buf;
	unsigned char* bytes = (unsigned char*)malloc(nbEle*sizeof(float));
	for(i=0;i<nbEle;i++)
	{
		buf.value = data[i];
		bytes[i*4+0] = buf.byte[0];
		bytes[i*4+1] = buf.byte[1];
		bytes[i*4+2] = buf.byte[2];
		bytes[i*4+3] = buf.byte[3];					
	}

	int byteLength = nbEle*sizeof(float);
	writeByteData(bytes, byteLength, tgtFilePath);
	free(bytes);
}

void writeDoubleData_inBytes(double *data, int nbEle, char* tgtFilePath)
{
	int i = 0, index = 0;
	ldouble buf;
	unsigned char* bytes = (unsigned char*)malloc(nbEle*sizeof(double));
	for(i=0;i<nbEle;i++)
	{
		index = i*8;
		buf.value = data[i];
		bytes[index+0] = buf.byte[0];
		bytes[index+1] = buf.byte[1];
		bytes[index+2] = buf.byte[2];
		bytes[index+3] = buf.byte[3];	
		bytes[index+4] = buf.byte[4];
		bytes[index+5] = buf.byte[5];
		bytes[index+6] = buf.byte[6];
		bytes[index+7] = buf.byte[7];							
	}

	int byteLength = nbEle*sizeof(double);
	writeByteData(bytes, byteLength, tgtFilePath);
	free(bytes);
}


void writeShortData(unsigned short *states, int stateLength, unsigned char *tgtFilePath)
{
	int byteLength = stateLength*2;
	unsigned char* bytes = (unsigned char*)malloc(byteLength*sizeof(char));
	convertShortArrayToBytes(states, stateLength, bytes);
	writeByteData(bytes, byteLength, tgtFilePath);
	free(bytes);
}

unsigned short* readShortData(unsigned char *srcFilePath, int *dataLength)
{
	int byteLength = 0;
	unsigned char * bytes = readByteData(srcFilePath, &byteLength);
	*dataLength = byteLength/2;
	unsigned short* states = convertByteDataToShortArray(bytes, byteLength);
	free(bytes);
	return states;
}
