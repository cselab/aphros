/**
 *  @file DynamicDoubleArray.h
 *  @author Sheng Di
 *  @date April, 2016
 *  @brief Header file for Dynamic Double Array.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _DynamicDoubleArray_H
#define _DynamicDoubleArray_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct DynamicDoubleArray
{	
	double* array;
	int size;
	int capacity;
} DynamicDoubleArray;

void new_DDA(DynamicDoubleArray **dda, int cap);
void convertDDAtoDoubles(DynamicDoubleArray *dba, double **data);
void free_DDA(DynamicDoubleArray *dda);
double getDDA_Data(DynamicDoubleArray *dda, int pos);
void addDDA_Data(DynamicDoubleArray *dda, int value);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _DynamicDoubleArray_H  ----- */
