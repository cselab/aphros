/**
 *  @file DynamicByteArray.h
 *  @author Sheng Di
 *  @date April, 2016
 *  @brief Header file for Dynamic Byte Array.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _DynamicByteArray_H
#define _DynamicByteArray_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct DynamicByteArray
{	
	unsigned char* array;
	int size;
	int capacity;
} DynamicByteArray;

void new_DBA(DynamicByteArray **dba, int cap);
void convertDBAtoBytes(DynamicByteArray *dba, unsigned char** bytes);
void free_DBA(DynamicByteArray *dba);
int getDBA_Data(DynamicByteArray *dba, int pos);
void addDBA_Data(DynamicByteArray *dba, unsigned char value);
void memcpyDBA_Data(DynamicByteArray *dba, unsigned char* data, int length);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _DynamicByteArray_H  ----- */
