/**
 *  @file io.h
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief Header file for the whole io interface.
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _IO_H
#define _IO_H

#include <stdio.h>
#include <stdio.h>

#ifdef _WIN32
#define PATH_SEPARATOR ';'
#else
#define PATH_SEPARATOR ':'
#endif

#ifdef __cplusplus
extern "C" {
#endif

int checkFileSize(char *srcFilePath);
unsigned char *readByteData(char *srcFilePath, int *byteLength);
double *readDoubleData_systemEndian(char *srcFilePath, int *nbEle);
float *readFloatData_systemEndian(char *srcFilePath, int *nbEle);
double *readDoubleData(char *srcFilePath, int *nbEle);
float *readFloatData(char *srcFilePath, int *nbEle);
void writeByteData(unsigned char *bytes, int outSize, char *tgtFilePath);
void writeDoubleData(double *data, int nbEle, char *tgtFilePath);
void writeFloatData(float *data, int nbEle, char *tgtFilePath);
void writeData(void *data, int dataType, int nbEle, char *tgtFilePath);
void writeFloatData_inBytes(float *data, int nbEle, char* tgtFilePath);
void writeDoubleData_inBytes(double *data, int nbEle, char* tgtFilePath);
void writeShortData(unsigned short *states, int stateLength, unsigned char *tgtFilePath);
unsigned short* readShortData(unsigned char *srcFilePath, int *dataLength);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _IO_H  ----- */
