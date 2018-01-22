/**
 *  @file test_decompress.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief This is an example of using Decompression interface.
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sz.h"
#include "rw.h"

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
    int nbEle, totalNbEle;
    char zipFilePath[640], outputFilePath[640];
    char *cfgFile;
    if(argc < 2)
    {
		printf("Test case: testfloat_decompress [configFile] [srcFilePath] [dimension sizes...]\n");
		printf("Example: testfloat_decompress sz.config testfloat_8_8_128.dat.sz 8 8 128\n");
		exit(0);
	}	
   
    cfgFile = argv[1];
    sprintf(zipFilePath, "%s", argv[2]);
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
    
    SZ_Init(cfgFile);
    
    if(r2==0)
	nbEle = r1;
    else if(r3==0)
	nbEle = r1*r2;
    else if(r4==0) 
    	nbEle = r1*r2*r3;
    else if(r5==0)
	nbEle = r1*r2*r3*r4;
    else
	nbEle = r1*r2*r3*r4*r5;

    sprintf(outputFilePath, "%s.out", zipFilePath);
    
    int byteLength;
    unsigned char *bytes = readByteData(zipFilePath, &byteLength);
  
    //printf("r1=%d,r2=%d,r3=%d,r4=%d,r5=%d\n", r1,r2,r3,r4,r5);
 
    cost_start();
    float *data = SZ_decompress(SZ_FLOAT, bytes, byteLength, r5, r4, r3, r2, r1);
    cost_end();
    //float data[r3][r2][r1];
    //nbEle = SZ_decompress_args(SZ_FLOAT, bytes, *byteLength, data, r5, r4, r3, r2, r1);
    
   // writeFloatData(data, nbEle, outputFilePath);
  
    free(bytes); 
    printf("timecost=%f\n",totalCost); 
    writeFloatData_inBytes(data, nbEle, outputFilePath);
    printf("done\n");
    
    SZ_Finalize();
    
    char oriFilePath[640];
    strncpy(oriFilePath, zipFilePath, (unsigned)strlen(zipFilePath)-3);
    oriFilePath[strlen(zipFilePath)-3] = '\0';
    float *ori_data = readFloatData(oriFilePath, &totalNbEle);
    int i = 0;
    float Max, Min, diffMax;
    Max = ori_data[0];
    Min = ori_data[0];
    diffMax = fabs(data[0] - ori_data[0]);
    int k = 0;
    double sum1 = 0, sum2 = 0;
    for (i = 0; i < nbEle; i++)
    {
        sum1 += ori_data[i];
	sum2 += data[i];
    }
    double mean1 = sum1/nbEle;
    double mean2 = sum2/nbEle;

    double sum3 = 0, sum4 = 0;
    double sum = 0, prodSum = 0;
    
    for (i = 0; i < nbEle; i++)
    {
        if (Max < ori_data[i]) Max = ori_data[i];
        if (Min > ori_data[i]) Min = ori_data[i];
        
        float err = fabs(data[i] - ori_data[i]);
        /*if(err>0.00001)
	{
		printf("%d %d: err=%.20G ori=%.20G dec=%.20G\n", k, i, err, ori_data[i], data[i]);
		break;
	}*/
	if (diffMax < err)
		diffMax = err;
        prodSum += (ori_data[i]-mean1)*(data[i]-mean2);
        sum3 += (ori_data[i] - mean1)*(ori_data[i]-mean1);
        sum4 += (data[i] - mean2)*(data[i]-mean2);
	sum += err*err;	
    }
    double std1 = sqrt(sum3/nbEle);
    double std2 = sqrt(sum4/nbEle);
    double ee = prodSum/nbEle;
    double acEff = ee/std1/std2;
 
    double mse = sum/nbEle;
    double range = Max - Min;
    double psnr = 20*log10(range)-10*log10(mse);
    double nrmse = sqrt(mse)/range;

    printf ("Min=%.20G, Max=%.20G, range=%.20G\n", Min, Max, range);
    printf ("Max absolute error = %.10f\n", diffMax);
    printf ("Max relative error = %f\n", diffMax/(Max-Min));
    printf ("PSNR = %f, NRMSE= %.20G\n", psnr,nrmse);
    printf ("acEff=%f\n", acEff);

    return 0;
}
