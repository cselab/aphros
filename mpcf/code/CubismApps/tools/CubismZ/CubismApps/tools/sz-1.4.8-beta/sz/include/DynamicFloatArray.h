/**
 *  @file DynamicFloatArray.h
 *  @author Sheng Di
 *  @date April, 2016
 *  @brief Header file for Dynamic Float Array.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _DynamicFloatArray_H
#define _DynamicFloatArray_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct DynamicFloatArray
{	
	float* array;
	int size;
	int capacity;
} DynamicFloatArray;

void new_DFA(DynamicFloatArray **dfa, int cap);
void convertDFAtoFloats(DynamicFloatArray *dfa, float **data);
void free_DFA(DynamicFloatArray *dfa);
float getDFA_Data(DynamicFloatArray *dfa, int pos);
void addDFA_Data(DynamicFloatArray *dfa, int value);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _DynamicFloatArray_H  ----- */
