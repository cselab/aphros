/**
 *  @file sz.c
 *  @author Sheng Di and Dingwen Tao
 *  @date Aug, 2016
 *  @brief SZ_Init, Compression and Decompression functions
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "sz.h"
#include "CompressElement.h"
#include "DynamicByteArray.h"
#include "DynamicIntArray.h"
#include "TightDataPointStorageD.h"
#include "TightDataPointStorageF.h"
#include "zlib.h"
#include "rw.h"
//#include "CurveFillingCompressStorage.h"

unsigned int maxRangeRadius = 32768;

int sysEndianType; //endian type of the system
int dataEndianType; //endian type of the data

char maxHeap[10];

long status;

int sol_ID;
int errorBoundMode; //ABS, REL, ABS_AND_REL, or ABS_OR_REL

int gzipMode; //four options: Z_NO_COMPRESSION, or Z_BEST_SPEED, Z_BEST_COMPRESSION, Z_DEFAULT_COMPRESSION

char *sz_cfgFile;

int offset;

double absErrBound;
double relBoundRatio;

int versionNumber[3];

int spaceFillingCurveTransform; //default is 0, or 1 set by sz.config
int reOrgSize; //the granularity of the reganization of the original data

int intvCapacity;
int intvRadius;

int layers = 1;
float predThreshold = 0.97;
int sampleDistance = 10;
char optQuantMode = 0; //opt Quantization (0: fixed ; 1: optimized)

int szMode = SZ_BEST_COMPRESSION;

SZ_VarSet* sz_varset;

sz_params *conf_params = NULL;

int SZ_Init(char *configFilePath)
{
	char str[512]="", str2[512]="", str3[512]="";
	sz_cfgFile = configFilePath;
	int loadFileResult = SZ_LoadConf();
	if(loadFileResult==0)
		exit(0);
	sz_varset = (SZ_VarSet*)malloc(sizeof(SZ_VarSet));
	sz_varset->header = (SZ_Variable*)malloc(sizeof(SZ_Variable));
	sz_varset->lastVar = sz_varset->header;
	sz_varset->count = 0;
	
	return 0;
}

void SZ_Reset()
{
	int i =0;
    if(pool==NULL)
    {
		pool = (struct node_t*)malloc(allNodes*2*sizeof(struct node_t));
		qqq = (node*)malloc(allNodes*2*sizeof(node));
		code = (unsigned long**)malloc(stateNum*sizeof(unsigned long*));//TODO
		cout = (unsigned char *)malloc(stateNum*sizeof(unsigned char));
	}
	
	memset(pool, 0, allNodes*2*sizeof(struct node_t));
	memset(qqq, 0, allNodes*2*sizeof(node));
    memset(code, 0, stateNum*sizeof(unsigned long*));
    memset(cout, 0, stateNum*sizeof(unsigned char));
	qq = qqq - 1;
	n_nodes = 0;
    n_inode = 0;
    qend = 1;
}

int SZ_Init_Params(sz_params *params)
{
	conf_params = params;
    int x = 1;
    char *y = (char*)&x;
    int endianType = BIG_ENDIAN_SYSTEM;
    if(*y==1) endianType = LITTLE_ENDIAN_SYSTEM;

    // set default values
    if(params->max_quant_intervals) 
    {
		maxRangeRadius = params->max_quant_intervals/2;
    	stateNum = maxRangeRadius*2;
		allNodes = maxRangeRadius*4;
		intvCapacity = maxRangeRadius*2;
		intvRadius = maxRangeRadius;		
    }
    dataEndianType    = endianType;
    sysEndianType    = endianType;
    sol_ID                    = SZ;
    offset                    = 0;
    gzipMode               = Z_BEST_SPEED;
    sampleDistance = 50;
    predThreshold = 0.97;
    errorBoundMode    = REL;
    absErrBound         = 0.000001;
    relBoundRatio         = 0.001;
    szMode = SZ_BEST_COMPRESSION;

    // set values from function arguments if avail.
    // [ENV]
    if(params->dataEndianType >= 0) dataEndianType    = params->dataEndianType;
    if(params->sysEndianType >= 0)    sysEndianType    = params->sysEndianType;
    if(params->sol_ID >= 0)                sol_ID                    = params->sol_ID;

    // [PARAMETER]
    if(sol_ID==SZ) {
        if(params->offset >= 0) offset = params->offset;

        /* gzipModes:
            Gzip_NO_COMPRESSION=0,
            Gzip_BEST_SPEED=1,
            Gzip_BEST_COMPRESSION=9,
            Gzip_DEFAULT_COMPRESSION=-1 */
        if(params->gzipMode >= -1) gzipMode = params->gzipMode;

		if(params->szMode >= 0) szMode = params->szMode;
        //if(params->maxSegmentNum >= 0) maxSegmentNum = params->maxSegmentNum;
        //if(params->spaceFillingCurveTransform >= 0) spaceFillingCurveTransform = params->spaceFillingCurveTransform;
        //if(params->reOrgSize >= 0) reOrgSize = params->reOrgSize;

        /* errBoundModes:
            ABS = 0
            REL = 1
            ABS_AND_REL =2
            ABS_OR_REL = 3 */
        if( params->errorBoundMode >= 0)  errorBoundMode =  params->errorBoundMode;

        if(params->absErrBound >= 0) absErrBound = params->absErrBound;
        if(params->relBoundRatio >= 0) relBoundRatio = params->relBoundRatio;
		if(params->quantization_intervals>0)
		{
			updateQuantizationInfo(params->quantization_intervals);
			optQuantMode = 0;
		}
		else
			optQuantMode = 1;
		if(params->layers >= 0)
			layers = params->layers;
		if(params->sampleDistance >= 0)
			sampleDistance = params->sampleDistance;
		if(params->predThreshold > 0)
			predThreshold = params->predThreshold;
    }

	versionNumber[0] = SZ_VER_MAJOR; //0
	versionNumber[1] = SZ_VER_MINOR; //5
	versionNumber[2] = SZ_VER_REVISION; //15

	if(params->quantization_intervals%2!=0)
	{
		printf("Error: quantization_intervals must be an even number!\n");
		return 0;
	}

    //initialization for Huffman encoding
	//SZ_Reset();
	
    return 1;
}

int computeDimension(int r5, int r4, int r3, int r2, int r1)
{
	int dimension;
	if(r1==0) 
	{
		dimension = 0;
	}
	else if(r2==0) 
	{
		dimension = 1;
	}
	else if(r3==0) 
	{
		dimension = 2;
	}
	else if(r4==0) 
	{
		dimension = 3;
	}
	else if(r5==0) 
	{
		dimension = 4;
	}
	else 
	{
		dimension = 5;
	}
	return dimension;	
}

int computeDataLength(int r5, int r4, int r3, int r2, int r1)
{
	int dataLength;
	if(r1==0) 
	{
		dataLength = 0;
	}
	else if(r2==0) 
	{
		dataLength = r1;
	}
	else if(r3==0) 
	{
		dataLength = r1*r2;
	}
	else if(r4==0) 
	{
		dataLength = r1*r2*r3;
	}
	else if(r5==0) 
	{
		dataLength = r1*r2*r3*r4;
	}
	else 
	{
		dataLength = r1*r2*r3*r4*r5;
	}
	return dataLength;
}

/*-------------------------------------------------------------------------*/
/**
    @brief      Perform Compression 
    @param      data           data to be compressed
    @param      outSize        the size (in bytes) after compression
    @param		r5,r4,r3,r2,r1	the sizes of each dimension (supporting only 5 dimensions at most in this version.
    @return     compressed data (in binary stream)

 **/
/*-------------------------------------------------------------------------*/
unsigned char* SZ_compress_args(int dataType, void *data, int *outSize, int errBoundMode, double absErrBound, 
double relBoundRatio, int r5, int r4, int r3, int r2, int r1)
{
	//TODO
	if(dataType==SZ_FLOAT)
	{
		unsigned char *newByteData;
		
		SZ_compress_args_float(&newByteData, (float *)data, r5, r4, r3, r2, r1, 
		outSize, errBoundMode, absErrBound, relBoundRatio);
		
		return newByteData;
	}
	else if(dataType==SZ_DOUBLE)
	{
		unsigned char *newByteData;
		SZ_compress_args_double(&newByteData, (double *)data, r5, r4, r3, r2, r1, 
		outSize, errBoundMode, absErrBound, relBoundRatio);
		
		return newByteData;
	}
	else
	{
		printf("Error: dataType can only be SZ_FLOAT or SZ_DOUBLE.\n");
		exit(1);
		return NULL;
	}
}

int SZ_compress_args2(int dataType, void *data, unsigned char* compressed_bytes, int *outSize, int errBoundMode, double absErrBound, double relBoundRatio, int r5, int r4, int r3, int r2, int r1)
{
	unsigned char* bytes = SZ_compress_args(dataType, data, outSize, errBoundMode, absErrBound, relBoundRatio, r5, r4, r3, r2, r1);
    memcpy(compressed_bytes, bytes, *outSize);
    free(bytes); 
	return 0;
}

unsigned char *SZ_compress(int dataType, void *data, int *outSize, int r5, int r4, int r3, int r2, int r1)
{	
	unsigned char *newByteData = SZ_compress_args(dataType, data, outSize, errorBoundMode, absErrBound, relBoundRatio, r5, r4, r3, r2, r1);
	return newByteData;
}

//////////////////
/*-------------------------------------------------------------------------*/
/**
    @brief      Perform Compression 
    @param      data           data to be compressed
    @param		reservedValue  the reserved value
    @param      outSize        the size (in bytes) after compression
    @param		r5,r4,r3,r2,r1	the sizes of each dimension (supporting only 5 dimensions at most in this version.
    @return     compressed data (in binary stream)

 **/
/*-------------------------------------------------------------------------*/
unsigned char *SZ_compress_rev_args(int dataType, void *data, void *reservedValue, int *outSize, int errBoundMode, double absErrBound, double relBoundRatio, int r5, int r4, int r3, int r2, int r1)
{
	int dataLength;
	unsigned char *newByteData;
	dataLength = computeDataLength(r5,r4,r3,r2,r1);
	//TODO
	printf("SZ compression with reserved data is TO BE DONE LATER.\n");
	exit(0);
	
	return newByteData;	
}

int SZ_compress_rev_args2(int dataType, void *data, void *reservedValue, unsigned char* compressed_bytes, int *outSize, int errBoundMode, double absErrBound, double relBoundRatio, int r5, int r4, int r3, int r2, int r1)
{
	unsigned char* bytes = SZ_compress_rev_args(dataType, data, reservedValue, outSize, errBoundMode, absErrBound, relBoundRatio, r5, r4, r3, r2, r1);
	memcpy(compressed_bytes, bytes, *outSize);
	free(bytes); //free(bytes) is removed , because of dump error at MIRA system (PPC architecture), fixed?
	return 0;
}

unsigned char *SZ_compress_rev(int dataType, void *data, void *reservedValue, int *outSize, int r5, int r4, int r3, int r2, int r1)
{
	int dataLength;

	unsigned char *newByteData;
	dataLength = computeDataLength(r5,r4,r3,r2,r1);	
	//TODO
	printf("SZ compression with reserved data is TO BE DONE LATER.\n");
	exit(0);
	
	return newByteData;
}

void *SZ_decompress(int dataType, unsigned char *bytes, int byteLength, int r5, int r4, int r3, int r2, int r1)
{
	if(dataType == SZ_FLOAT)
	{
		float *newFloatData;
		SZ_decompress_args_float(&newFloatData, r5, r4, r3, r2, r1, bytes, byteLength);
		return newFloatData;
	}
	else if(dataType == SZ_DOUBLE)
	{
		double *newDoubleData;
		SZ_decompress_args_double(&newDoubleData, r5, r4, r3, r2, r1, bytes, byteLength);
		return newDoubleData;
	}
	else
	{
		printf("Error: data type cannot be the types other than SZ_FLOAT or SZ_DOUBLE\n");
		exit(0);		
	}
}

int SZ_decompress_args(int dataType, unsigned char *bytes, int byteLength, void* decompressed_array, int r5, int r4, int r3, int r2, int r1)
{
	int i;
	int nbEle = computeDataLength(r5,r4,r3,r2,r1);
	if(dataType == SZ_FLOAT)
	{
		float* data = (float *)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		float* data_array = (float *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(float));
		//for(i=0;i<nbEle;i++)
		//	data_array[i] = data[i];	
		free(data); //this free operation seems to not work with BlueG/Q system.
	}
	else if(dataType == SZ_DOUBLE)
	{
		double* data = (double *)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		double* data_array = (double *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(double));
		//for(i=0;i<nbEle;i++)
		//	data_array[i] = data[i];
		free(data); //this free operation seems to not work with BlueG/Q system.
	}
	else
	{
		printf("Error: data type cannot be the types other than SZ_FLOAT or SZ_DOUBLE\n");
		exit(0);				
	}

	return nbEle;
}

/*-----------------------------------batch data compression--------------------------------------*/

void filloutDimArray(int* dim, int r5, int r4, int r3, int r2, int r1)
{
	if(r2==0)
		dim[0] = r1;
	else if(r3==0)
	{
		dim[0] = r2;
		dim[1] = r1;
	}
	else if(r4==0)
	{
		dim[0] = r3;
		dim[1] = r2;
		dim[2] = r1;
	}
	else if(r5==0)
	{
		dim[0] = r4;
		dim[1] = r3;
		dim[2] = r2;
		dim[3] = r1;
	}
	else
	{
		dim[0] = r5;
		dim[1] = r4;
		dim[2] = r3;
		dim[3] = r2;
		dim[4] = r1;		
	}
}

unsigned char* SZ_batch_compress(int *outSize)
{
	int dataLength;
	DynamicByteArray* dba; 
	new_DBA(&dba, 32768);
	
	//number of variables
	int varCount = sz_varset->count;
	unsigned char countBufBytes[4];
	intToBytes_bigEndian(countBufBytes, varCount);
	//add only the lats two bytes, i.e., the maximum # variables is supposed to be less than 32768
	addDBA_Data(dba, countBufBytes[2]);
	addDBA_Data(dba, countBufBytes[3]);
	
	int i, j, k = 0;
	SZ_Variable* p = sz_varset->header->next;
	while(p!=NULL)
	{
		if(p->dataType==SZ_FLOAT)
		{
			unsigned char *newByteData;
			int outSize;
			SZ_compress_args_float_wRngeNoGzip(&newByteData, (float *)p->data, 
			p->r5, p->r4, p->r3, p->r2, p->r1, 
			&(p->compressedSize), p->errBoundMode, p->absErrBound, p->relBoundRatio);
			
			p->compressedBytes = newByteData;
		}
		else if(p->dataType==SZ_DOUBLE)
		{
			unsigned char *newByteData;
			int outSize;
			SZ_compress_args_double_wRngeNoGzip(&newByteData, (double *)p->data, 
			p->r5, p->r4, p->r3, p->r2, p->r1, 
			&(p->compressedSize), p->errBoundMode, p->absErrBound, p->relBoundRatio);
			
			p->compressedBytes = newByteData;
		}
		
		//TODO: copy metadata information (totally 1+4*x+4+(1+y) bytes)
		//byte format of variable storage: meta 1 bytes (000000xy) x=0 SZ/1 HZ, y=0 float/1 double ; 
		//4*x: x integers indicating the dimension sizes; 
		//compressed length (int) 4 bytes;
		//1+y: 1 means the length of varname, y bytes represents the varName string. 
		int meta = 0;
		if(p->dataType == SZ_FLOAT)
			meta = 2; //10
		else if(p->dataType == SZ_DOUBLE)
			meta = 3; //11
		
		//keep dimension information
		int dimNum = computeDimension(p->r5, p->r4, p->r3, p->r2, p->r1);
		int dimSize[dimNum];
		memset(dimSize, 0, dimNum*sizeof(int));
		meta = meta | dimNum << 2; //---aaabc: aaa indicates dim, b indicates HZ, c indicates dataType
		
		addDBA_Data(dba, (unsigned char)meta);
		
		filloutDimArray(dimSize, p->r5, p->r4, p->r3, p->r2, p->r1);
		
		for(j=0;j<dimNum;j++)
		{
			intToBytes_bigEndian(countBufBytes, dimSize[j]);
			for(i = 0;i<4;i++)
				addDBA_Data(dba, countBufBytes[i]);
		}
			 
		//Keep compressed size information	 
		intToBytes_bigEndian(countBufBytes, p->compressedSize);
		
		for(i = 0;i<4;i++)
			addDBA_Data(dba, countBufBytes[i]);			 
			 
		//Keep varName information	 
		unsigned char varNameLength = (unsigned char)strlen(p->varName);
		addDBA_Data(dba, varNameLength);
		memcpyDBA_Data(dba, (unsigned char*)p->varName, varNameLength);
			 
		//copy the compressed stream
		memcpyDBA_Data(dba, p->compressedBytes, p->compressedSize);
		free(p->compressedBytes);
		p->compressedBytes = NULL;
			 
		p = p->next;
	}
	
	unsigned char* tmpFinalCompressedBytes;
	convertDBAtoBytes(dba, &tmpFinalCompressedBytes);
	
	unsigned char* tmpCompressedBytes2;
	int tmpGzipSize = (int)zlib_compress2(tmpFinalCompressedBytes, dba->size, &tmpCompressedBytes2, gzipMode);
	free(tmpFinalCompressedBytes);
	
	unsigned char* finalCompressedBytes = (unsigned char*) malloc(sizeof(unsigned char)*(4+tmpGzipSize));
	
	intToBytes_bigEndian(countBufBytes, dba->size);
	
	memcpy(finalCompressedBytes, countBufBytes, 4);
	
	memcpy(&(finalCompressedBytes[4]), tmpCompressedBytes2, tmpGzipSize);
	free(tmpCompressedBytes2);
	
	*outSize = 4+tmpGzipSize;
	free(dba);
	return finalCompressedBytes;
}

SZ_VarSet* SZ_batch_decompress(unsigned char* compressedStream, int compressedLength)
{
	int i, j, k = 0;
	unsigned char intByteBuf[4];
	
	//get target decompression size for Gzip (zlib)
	intByteBuf[0] = compressedStream[0];
	intByteBuf[1] = compressedStream[1];
	intByteBuf[2] = compressedStream[2];
	intByteBuf[3] = compressedStream[3];
	
	int targetUncompressSize = bytesToInt_bigEndian(intByteBuf);
	
	//Gzip decompression
	unsigned char* gzipDecpressBytes;
	int gzipDecpressSize = zlib_uncompress2(&(compressedStream[4]), (unsigned long)compressedLength, &gzipDecpressBytes, (unsigned long)targetUncompressSize);
	
	if(gzipDecpressSize!=targetUncompressSize)
	{
		printf("Error: incorrect decompression in zlib_uncompress3: gzipDecpressSize!=targetUncompressSize\n");
		exit(0);
	}
	
	//Start analyzing the byte stream for further decompression	
	intByteBuf[0] = 0;
	intByteBuf[1] = 0; 
	intByteBuf[2] = gzipDecpressBytes[k++];
	intByteBuf[3] = gzipDecpressBytes[k++];
	
	int varCount = bytesToInt_bigEndian(intByteBuf);	
		
	int varNum = sz_varset->count;
	int dataLength, cpressedLength;
	
	if(varNum==0)
	{
		SZ_Variable* lastVar = sz_varset->header; 
		if(lastVar->next!=NULL)
		{
			printf("Error: sz_varset.count is inconsistent with the number of variables in sz_varset->header.\n");
			exit(0);
		}
		for(i=0;i<varCount;i++)
		{
			int type = (int)gzipDecpressBytes[k++];
			int dataType;
			int tmpType = type & 0x03;
			if(tmpType==2) //FLOAT
				dataType = SZ_FLOAT;
			else if(tmpType==3)//DOUBLE
				dataType = SZ_DOUBLE;
			else
			{
				printf("Error: Wrong value of decompressed type. \n Please check the correctness of the decompressed data.\n");
				exit(0);
			}
			
			//get # dimensions and the size of each dimension
			int dimNum = (type & 0x1C) >> 2; //compute dimension
			int start_dim = 5 - dimNum;
			int dimSize[5];
			memset(dimSize, 0, 5*sizeof(int));
			
			for(j=0;j<dimNum;j++)
			{
				memcpy(intByteBuf, &(gzipDecpressBytes[k]), 4);
				k+=4;
				dimSize[start_dim+j] = bytesToInt_bigEndian(intByteBuf);
			}	
			
			//get compressed length
			memcpy(intByteBuf, &(gzipDecpressBytes[k]), 4);
			k+=4;
			cpressedLength = bytesToInt_bigEndian(intByteBuf);	
					
			//Keep varName information	 
			int varNameLength = gzipDecpressBytes[k++];
			char* varNameString = (char*)malloc(sizeof(char)*(varNameLength+1));	
			memcpy(varNameString, &(gzipDecpressBytes[k]), varNameLength);
			k+=varNameLength;
			varNameString[varNameLength] = '\0';
		
			//TODO: convert szTmpBytes to data array.
			dataLength = computeDataLength(dimSize[0], dimSize[1], dimSize[2], dimSize[3], dimSize[4]);
			
			if(dataType==SZ_FLOAT)
			{
				float* newData;
				TightDataPointStorageF* tdps;
				new_TightDataPointStorageF_fromFlatBytes(&tdps, &(gzipDecpressBytes[k]), cpressedLength);

				if (dataLength == dimSize[4])
					getSnapshotData_float_1D(&newData,dimSize[4],tdps);
				else
				if (dataLength == dimSize[3]*dimSize[4])
					getSnapshotData_float_2D(&newData,dimSize[3],dimSize[4],tdps);
				else
				if (dataLength == dimSize[2]*dimSize[3]*dimSize[4])
					getSnapshotData_float_3D(&newData,dimSize[2], dimSize[3], dimSize[4],tdps);
				else
				if (dataLength == dimSize[1]*dimSize[2]*dimSize[3]*dimSize[4])
					getSnapshotData_float_3D(&newData,dimSize[1]*dimSize[2], dimSize[3], dimSize[4],tdps);
				else
				{
					printf("Current version doesn't support 5 dimensions.\n");
					exit(0);
				}
				
				SZ_batchAddVar(varNameString, dataType, newData, dimSize[0], dimSize[1], dimSize[2], dimSize[3], dimSize[4], 0, 0, 0);
				free_TightDataPointStorageF(tdps);			
			}	
			else if(dataType==SZ_DOUBLE)
			{
				double* newData;
				TightDataPointStorageD* tdps;
				new_TightDataPointStorageD_fromFlatBytes(&tdps, &(gzipDecpressBytes[k]), cpressedLength);

				if (dataLength == dimSize[4])
					getSnapshotData_double_1D(&newData,dimSize[4],tdps);
				else
				if (dataLength == dimSize[3]*dimSize[4])
					getSnapshotData_double_2D(&newData,dimSize[3],dimSize[4],tdps);
				else
				if (dataLength == dimSize[2]*dimSize[3]*dimSize[4])
					getSnapshotData_double_3D(&newData,dimSize[2], dimSize[3], dimSize[4],tdps);
				else
				if (dataLength == dimSize[1]*dimSize[2]*dimSize[3]*dimSize[4])
					getSnapshotData_double_3D(&newData,dimSize[1]*dimSize[2], dimSize[3], dimSize[4],tdps);
				else
				{
					printf("Current version doesn't support 5 dimensions.\n");
					exit(0);
				}
			
				SZ_batchAddVar(varNameString, dataType, newData, dimSize[0], dimSize[1], dimSize[2], dimSize[3], dimSize[4], 0, 0, 0);
				free_TightDataPointStorageD(tdps);					
			}
			else
			{
				printf("Error: wrong dataType in the batch decompression.\n");
				exit(0);
			}
	
			k+=cpressedLength;
		}
	}
	else if(varNum!=varCount)
	{
		printf("Error: Inconsistency of the # variables between sz_varset and decompressed stream.\n");
		printf("Note sz_varset.count should be either 0 or the correct number of variables stored in the decompressed stream.\n");
		printf("Currently, sz_varset.count=%d, expected number of variables = %d\n", varNum, varCount);
		exit(0);
	}
	else //if(varNum>0 && varNum==varCount)
	{
		SZ_Variable* p = sz_varset->header; 
		for(i=0;i<varCount;i++)
		{
			int type = (int)gzipDecpressBytes[k++];
			int dataType;
			int tmpType = type & 0x03;
			if(tmpType==2) //FLOAT
				dataType = SZ_FLOAT;
			else if(tmpType==3)//DOUBLE
				dataType = SZ_DOUBLE;
			else
			{
				printf("Error: Wrong value of decompressed type. \n Please check the correctness of the decompressed data.\n");
				exit(0);
			}
			
			//get # dimensions and the size of each dimension
			int dimNum = (type & 0x1C) >> 2; //compute dimension
			int start_dim = 5 - dimNum;
			int dimSize[5];
			memset(dimSize, 0, 5*sizeof(int));
			
			for(j=0;j<dimNum;j++)
			{
				memcpy(intByteBuf, &(gzipDecpressBytes[k]), 4);
				k+=4;
				dimSize[start_dim+j] = bytesToInt_bigEndian(intByteBuf);
			}
			
			//get compressed length
			memcpy(intByteBuf, &(gzipDecpressBytes[k]), 4);
			k+=4;
			cpressedLength = bytesToInt_bigEndian(intByteBuf);	
			
			//Keep varName information	 
			int varNameLength = gzipDecpressBytes[k++];
			char* varNameString = (char*)malloc(sizeof(char)*(varNameLength+1));	
			memcpy(varNameString, &(gzipDecpressBytes[k]), varNameLength);
			k+=varNameLength;
			varNameString[varNameLength] = '\0';
			
			SZ_Variable* curVar = p->next;
			
			int checkVarName = strcmp(varNameString, curVar->varName);
			if(checkVarName!=0)
			{
				printf("Error: the varNames in the compressed stream are inconsistent with those in the sz_varset\n");
				exit(100);
			}
						
			//TODO: convert szTmpBytes to data array.
			dataLength = computeDataLength(dimSize[0], dimSize[1], dimSize[2], dimSize[3], dimSize[4]);
			
			if(dataType==SZ_FLOAT)
			{
				float* newData;
				TightDataPointStorageF* tdps;
				new_TightDataPointStorageF_fromFlatBytes(&tdps, &(gzipDecpressBytes[k]), cpressedLength);					

				if (dataLength == dimSize[4])
					getSnapshotData_float_1D(&newData,dimSize[4],tdps);
				else
				if (dataLength == dimSize[3]*dimSize[4])
					getSnapshotData_float_2D(&newData,dimSize[3],dimSize[4],tdps);
				else
				if (dataLength == dimSize[2]*dimSize[3]*dimSize[4])
					getSnapshotData_float_3D(&newData,dimSize[2], dimSize[3], dimSize[4],tdps);
				else
				if (dataLength == dimSize[1]*dimSize[2]*dimSize[3]*dimSize[4])
					getSnapshotData_float_3D(&newData,dimSize[1]*dimSize[2], dimSize[3], dimSize[4],tdps);
				else
				{
					printf("Error: doesn't support 5 dimensions yet.\n");
					exit(0);
				}
				
				free_TightDataPointStorageF(tdps);	
				curVar->data = newData;
			}	
			else if(dataType==SZ_DOUBLE)
			{
				double* newData;
				TightDataPointStorageD* tdps;
				new_TightDataPointStorageD_fromFlatBytes(&tdps, &(gzipDecpressBytes[k]), cpressedLength);

				if (dataLength == dimSize[4])
					getSnapshotData_double_1D(&newData,dimSize[4],tdps);
				else
				if (dataLength == dimSize[3]*dimSize[4])
					getSnapshotData_double_2D(&newData,dimSize[3],dimSize[4],tdps);
				else
				if (dataLength == dimSize[2]*dimSize[3]*dimSize[4])
					getSnapshotData_double_3D(&newData,dimSize[2], dimSize[3], dimSize[4],tdps);
				else
				if (dataLength == dimSize[1]*dimSize[2]*dimSize[3]*dimSize[4])
					getSnapshotData_double_3D(&newData,dimSize[1]*dimSize[2], dimSize[3], dimSize[4],tdps);
				else
				{
					printf("Error: doesn't support 5 dimensions yet.\n");
					exit(0);
				}
				SZ_batchAddVar(varNameString, dataType, newData, dimSize[0], dimSize[1], dimSize[2], dimSize[3], dimSize[4], 0, 0, 0);
				free_TightDataPointStorageD(tdps);	
			}
			else
			{
				printf("Error: wrong dataType in the batch decompression.\n");
				exit(0);
			}		
			
			free(varNameString);
			k+=cpressedLength;
			p = p->next;
		}
	}
	free(gzipDecpressBytes);
	return sz_varset;
}

/**
 * @deprecated
 * @return: the length of the coefficient array.
 * */
int getPredictionCoefficients(int layers, int dimension, int **coeff_array)
{
	int size = 0;
	switch(dimension)
	{
		case 1:
			switch(layers)
			{
				case 1:
					*coeff_array = (int*)malloc(sizeof(int));
					(*coeff_array)[0] = 1;
					size = 1;
					break;
				case 2:
					*coeff_array = (int*)malloc(2*sizeof(int));
					(*coeff_array)[0] = 2;
					(*coeff_array)[1] = -1;
					size = 2;
					break;
				case 3:
					*coeff_array = (int*)malloc(3*sizeof(int));
					(*coeff_array)[0] = 3;
					(*coeff_array)[1] = -3;
					(*coeff_array)[2] = 1;
					break;
			}	
			break;
		case 2:
			switch(layers)
			{
				case 1:
				
					break;
				case 2:
				
					break;
				case 3:
				
					break;
			}				
			break;
		case 3:
			switch(layers)
			{
				case 1:
				
					break;
				case 2:
				
					break;
				case 3:
				
					break;
			}			
			break;
		default:
			printf("Error: dimension must be no greater than 3 in the current version.\n");
			exit(1);
	}
	return size;
}

void SZ_Finalize()
{
	int i ;
	free(pool);
	pool = NULL;
	free(qqq);
	qqq = NULL;
	for(i=0;i<stateNum;i++)
	{
		if(code[i]!=NULL)
			free(code[i]);
	}	
	free(code);
	code = NULL;
	free(cout);
	cout = NULL;
	
}
