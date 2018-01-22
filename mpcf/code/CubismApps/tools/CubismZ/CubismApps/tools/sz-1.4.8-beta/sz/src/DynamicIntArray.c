/**
 *  @file DynamicIntArray.c
 *  @author Sheng Di
 *  @date May, 2016
 *  @brief Dynamic Int Array
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include "DynamicIntArray.h"

void new_DIA(DynamicIntArray **dia, int cap) {
		*dia = (DynamicIntArray *)malloc(sizeof(DynamicIntArray));
        (*dia)->size = 0;
        (*dia)->capacity = cap;
        (*dia)->array = (unsigned char*)malloc(sizeof(unsigned char)*cap);
    }

void convertDIAtoInts(DynamicIntArray *dia, unsigned char **data)
{
	int size = dia->size;
	if(size>0)
		*data = (unsigned char*)malloc(size * sizeof(char));
	else
		*data = NULL;
	memcpy(*data, dia->array, size*sizeof(unsigned char));	
}

void free_DIA(DynamicIntArray *dia)
{
	free(dia->array);
	free(dia);
}

int getDIA_Data(DynamicIntArray *dia, int pos)
{
	if(pos>=dia->size)
	{
		printf("Error: wrong position of DIA.\n");
		exit(0);
	}
	return dia->array[pos];
}

void addDIA_Data(DynamicIntArray *dia, int value)
{
	if(dia->size==dia->capacity)
	{
		dia->capacity = dia->capacity << 1;
		dia->array = (unsigned char *)realloc(dia->array, dia->capacity*sizeof(unsigned char));
	}
	dia->array[dia->size] = (unsigned char)value;
	dia->size ++;
}
