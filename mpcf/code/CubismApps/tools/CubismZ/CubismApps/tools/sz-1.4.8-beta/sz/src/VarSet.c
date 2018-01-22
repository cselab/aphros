/**
 *  @file Variable.c
 *  @author Sheng Di
 *  @date July, 2016
 *  @brief TypeManager is used to manage the type array: parsing of the bytes and other types in between.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "VarSet.h"
#include "sz.h"

void free_Variable_keepDecompressedData(SZ_Variable* v)
{
	if(v->data!=NULL)
		free(v->compressedBytes);
	
	free(v);
}

void free_Variable_keepCompressedBytes(SZ_Variable* v)
{
	if(v->data!=NULL)
		free(v->data);
	
	free(v);
}

void free_Variable_all(SZ_Variable* v)
{
	if(v->data!=NULL)
		free(v->data);
	if(v->compressedBytes!=NULL)
		free(v->compressedBytes);
	
	free(v);
}

void SZ_batchAddVar(char* varName, int dataType, void* data, 
			int errBoundMode, double absErrBound, double relBoundRatio,
			int r5, int r4, int r3, int r2, int r1)
{	
	SZ_Variable* var = (SZ_Variable*)malloc(sizeof(SZ_Variable));
	
	var->varName = varName;
	var->dataType = dataType;
	var->r5 = r5;
	var->r4 = r4;
	var->r3 = r3;
	var->r2 = r2;
	var->r1 = r1;
	var->errBoundMode = errBoundMode;
	var->absErrBound = absErrBound;
	var->relBoundRatio = relBoundRatio;
	var->data = data;
	var->compressedBytes = NULL;
	var->next = NULL;
	
	sz_varset->count ++;
	sz_varset->lastVar->next = var;
	sz_varset->lastVar = var;
}

int SZ_batchDelVar(char* varName)
{
	int state = SZ_batchDelVar_vset(sz_varset, varName);
	return state;
}

int SZ_batchDelVar_vset(SZ_VarSet* vset, char* varName)
{
	int delSuccess = 1;
	SZ_Variable* p = vset->header;
	SZ_Variable* q = p->next;
	while(q != NULL)
	{
		int cmpResult = strcmp(q->varName, varName);
		if(cmpResult==0)
		{
			p->next = q->next;
			free_Variable_all(q);
			vset->count --;
			delSuccess = 0;
			break;
		}
		p = p->next;
		q = q->next;	
	}
	
	return delSuccess;
}

SZ_Variable* SZ_searchVar(char* varName)
{
	SZ_Variable* p = sz_varset->header->next;
	while(p!=NULL)
	{
		int checkName = strcmp(p->varName, varName);
		if(checkName==0)
			return p;
		p = p->next;
	}	
	return NULL;
}

void* SZ_getVarData(char* varName, int *r5, int *r4, int *r3, int *r2, int *r1)
{
	SZ_Variable* v = SZ_searchVar(varName);
	*r5 = v->r5;
	*r4 = v->r4;
	*r3 = v->r3;
	*r2 = v->r2;
	*r1 = v->r1;
	return (void*)v->data;
}

void free_VarSet(void)
{
	free_VarSet_vset(sz_varset);
}

//free_VarSet will completely destroy the SZ_VarSet, so don't do it until you really don't need to any more!
void free_VarSet_vset(SZ_VarSet *vset)
{
	SZ_Variable *p = vset->header;
	while(p->next!=NULL)
	{
		SZ_Variable *q = p->next;
		p->next = q->next;
		free_Variable_all(q);
	}
	free(vset);
}
