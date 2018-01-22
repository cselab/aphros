/**
 *  @file testdouble_batch_compress.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief This is an example of using compression interface
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include "sz.h"
#include "rw.h"
#include "VarSet.h"

struct timeval startTime;
struct timeval endTime;  /* Start and end times */
struct timeval costStart; /*only used for recording the cost*/
double totalCost = 0;


void cost_start()
{
        gettimeofday(&costStart, NULL);
}

void cost_end()
{
        double elapsed;
        struct timeval costEnd;
        gettimeofday(&costEnd, NULL);
        elapsed = ((costEnd.tv_sec*1000000+costEnd.tv_usec)-(costStart.tv_sec*1000000+costStart.tv_usec))/1000000.0;
        totalCost += elapsed;
}


int main(int argc, char * argv[])
{
    int r5=0,r4=0,r3=0,r2=0,r1=0;
    char outDir[640], oriFilePath[640], outputFilePath[640];
    char *cfgFile;
    
    if(argc < 3)
    {
	printf("Test case: testdouble_batch_compress [config_file] [srcFilePath] [dimension sizes...]\n");
	printf("Example: testdouble_batch_compress sz.config testdouble_8_8_128.dat 8 8 128\n");
	exit(0);
    }
   
    cfgFile=argv[1];
    sprintf(oriFilePath, "%s", argv[2]);
    if(argc>=4)
	r1 = atoi(argv[3]); //8
    if(argc>=5)
	r2 = atoi(argv[4]); //8
    if(argc>=6)
	r3 = atoi(argv[5]); //128
    if(argc>=7)
        r4 = atoi(argv[6]);
    if(argc>=8)
        r5 = atoi(argv[7]);
   
    printf("cfgFile=%s\n", cfgFile); 
    SZ_Init(cfgFile);
    
    sprintf(outputFilePath, "%s.bsz", oriFilePath);
   
    int nbEle;
    double *data = readDoubleData(oriFilePath, &nbEle);
  
    double *data2 = readDoubleData(oriFilePath, &nbEle);
    
    SZ_batchAddVar("data", SZ_DOUBLE, data, ABS, 0.000001, 0.000001, r5,r4,r3,r2,r1);
    SZ_batchAddVar("data2", SZ_DOUBLE, data2, ABS, 0.000001, 0.000001, r5,r4,r3,r2,r1);

    int outSize; 
    cost_start();
    unsigned char *bytes = SZ_batch_compress(&outSize);
    cost_end(); 
    printf("timecost=%f\n",totalCost); 
    writeByteData(bytes, outSize, outputFilePath);
   
    //decompression
    int compressedLength = outSize;
    SZ_batch_decompress(bytes, compressedLength);

    //write the decompressed data to disk 
    int dataLength, varIndex = 0;
    SZ_Variable* p = sz_varset->header->next;
    while(p!=NULL)
    {
        sprintf(outputFilePath, "%s-batch%d.out", oriFilePath, varIndex++);
        dataLength = computeDataLength(p->r5, p->r4, p->r3, p->r2, p->r1);
	writeDoubleData(p->data, dataLength, outputFilePath);
        p = p->next;
    }

    printf("done\n");
    
    free(sz_varset);
    SZ_Finalize();
    
    return 0;
}
