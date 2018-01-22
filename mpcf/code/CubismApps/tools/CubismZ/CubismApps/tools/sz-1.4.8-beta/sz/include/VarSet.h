/**
 *  @file VarSet.h
 *  @author Sheng Di
 *  @date July, 2016
 *  @brief Header file for the Variable.c.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _VarSet_H
#define _VarSet_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SZ_Variable
{
	char* varName;
	char compressType; //102 means HZ; 101 means SZ 
	int dataType; //SZ_FLOAT or SZ_DOUBLE
	int r5;
	int r4;
	int r3;
	int r2;
	int r1;
	int errBoundMode;
	double absErrBound;
	double relBoundRatio;
	void* data;
	unsigned char* compressedBytes;
	int compressedSize;
	struct SZ_Variable* next;
} SZ_Variable;

typedef struct SZ_VarSet
{
	int count;
	struct SZ_Variable *header;
	struct SZ_Variable *lastVar;
} SZ_VarSet;

void free_Variable_keepDecompressedData(SZ_Variable* v);
void free_Variable_keepCompressedBytes(SZ_Variable* v);
void free_Variable_all(SZ_Variable* v);
void SZ_batchAddVar(char* varName, int dataType, void* data, 
			int errBoundMode, double absErrBound, double relBoundRatio,
			int r5, int r4, int r3, int r2, int r1);
int SZ_batchDelVar_vset(SZ_VarSet* vset, char* varName);
int SZ_batchDelVar(char* varName);

SZ_Variable* SZ_searchVar(char* varName);
void* SZ_getVarData(char* varName, int *r5, int *r4, int *r3, int *r2, int *r1);

void free_VarSet_vset(SZ_VarSet *vset);
void free_VarSet(void);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _VarSet_H  ----- */
