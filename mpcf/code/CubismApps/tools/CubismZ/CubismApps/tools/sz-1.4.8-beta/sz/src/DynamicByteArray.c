/**
 *  @file DynamicByteArray.c
 *  @author Sheng Di
 *  @date May, 2016
 *  @brief Dynamic Byte Array
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include "DynamicByteArray.h"

void new_DBA(DynamicByteArray **dba, int cap) {
		*dba = (DynamicByteArray *)malloc(sizeof(DynamicByteArray));
        (*dba)->size = 0;
        (*dba)->capacity = cap;
        (*dba)->array = (unsigned char*)malloc(sizeof(unsigned char)*cap);
    }

void convertDBAtoBytes(DynamicByteArray *dba, unsigned char** bytes)
{
	int size = dba->size;
	if(size>0)
		*bytes = (unsigned char*)malloc(size * sizeof(unsigned char));
	else
		*bytes = NULL;
	memcpy(*bytes, dba->array, size*sizeof(unsigned char));	
}

void free_DBA(DynamicByteArray *dba)
{
	free(dba->array);
	free(dba);
}

int getDBA_Data(DynamicByteArray *dba, int pos)
{
	if(pos>=dba->size)
	{
		printf("Error: wrong position of DBA.\n");
		exit(0);
	}
	return dba->array[pos];
}

void addDBA_Data(DynamicByteArray *dba, unsigned char value)
{
	if(dba->size==dba->capacity)
	{
		dba->capacity = dba->capacity << 1;
		dba->array = (unsigned char *)realloc(dba->array, dba->capacity*sizeof(unsigned char));
	}
	dba->array[dba->size] = value;
	dba->size ++;
}

void memcpyDBA_Data(DynamicByteArray *dba, unsigned char* data, int length)
{
	if(dba->size + length > dba->capacity)
	{
		dba->capacity = dba->size + length;
		dba->array = (unsigned char *)realloc(dba->array, dba->capacity*sizeof(unsigned char));
	}
	memcpy(&(dba->array[dba->size]), data, length);
	dba->size += length;
}
