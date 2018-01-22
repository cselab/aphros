/**
 *  @file   conf.c
 *  @author Sheng Di (sdi1@anl.gov or disheng222@gmail.com)
 *  @date   2015.
 *  @brief  Configuration loading functions for the SZ library.
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "string.h"
#include "sz.h"
#include "iniparser.h"
#include "Huffman.h"

/*-------------------------------------------------------------------------*/
/**
    @brief      It reads the configuration given in the configuration file.
    @return     integer         1 if successfull.

    This function reads the configuration given in the SZ configuration
    file and sets other required parameters.

 **/
 
/*struct node_t *pool;
node *qqq;
node *qq;
int n_nodes = 0, qend;
unsigned long **code;
unsigned char *cout;
int n_inode;*/ 
 
unsigned int roundUpToPowerOf2(unsigned int base)
{
  base -= 1;

  base = base | (base >> 1);
  base = base | (base >> 2);
  base = base | (base >> 4);
  base = base | (base >> 8);
  base = base | (base >> 16);

  return base + 1;
} 
 
void updateQuantizationInfo(int quant_intervals)
{
	allNodes = 2*quant_intervals;
	stateNum = quant_intervals;
	intvCapacity = quant_intervals;
	intvRadius = quant_intervals/2;
} 
 
void clearHuffmanMem()
{
   	memset(pool, 0, allNodes*sizeof(struct node_t));
	memset(code, 0, stateNum*sizeof(long*)); //original:128; we just used '0'-'7', so max ascii is 55.
	memset(cout, 0, stateNum);
    n_nodes = 0;
    n_inode = 0;
    qend = 1;	
} 
 
/*-------------------------------------------------------------------------*/
int SZ_ReadConf() {
    // Check access to SZ configuration file and load dictionary
    //record the setting in conf_params
    conf_params = (sz_params*)malloc(sizeof(sz_params));    
    
    int x = 1;
    char sol_name[256];
    char *modeBuf;
    char *errBoundMode;
    char *endianTypeString;
    dictionary *ini;
    char *par;
    printf("[SZ] Reading SZ configuration file (%s) ...\n", sz_cfgFile);
    if (access(sz_cfgFile, F_OK) != 0)
    {
        printf("[SZ] Configuration file NOT accessible.\n");
        return 1;
    }
    ini = iniparser_load(sz_cfgFile);
    if (ini == NULL)
    {
        printf("[SZ] Iniparser failed to parse the conf. file.\n");
        return 1;
    }

	endianTypeString = iniparser_getstring(ini, "ENV:dataEndianType", NULL);
	if(strcmp(endianTypeString, "LITTLE_ENDIAN_DATA")==0)
		dataEndianType = LITTLE_ENDIAN_DATA;
	else if(strcmp(endianTypeString, "BIG_ENDIAN_DATA")==0)
		dataEndianType = BIG_ENDIAN_DATA;
	else
	{
		printf("Error: Wrong dataEndianType: please set it correctly in sz.config.\n");
		iniparser_freedict(ini);
		exit(1);
	}

	conf_params->dataEndianType = dataEndianType;

	char *y = (char*)&x;
	
	if(*y==1)
		sysEndianType = LITTLE_ENDIAN_SYSTEM;
	else //=0
		sysEndianType = BIG_ENDIAN_SYSTEM;
	conf_params->sysEndianType = sysEndianType;

	// Reading/setting detection parameters
	
	par = iniparser_getstring(ini, "ENV:sol_name", NULL);
	snprintf(sol_name, 256, "%s", par);
	
    if(strcmp(sol_name, "SZ")==0)
		sol_ID = SZ;
	else
	{
		printf("[SZ] Error: wrong solution name (please check sz.config file)\n");
		iniparser_freedict(ini);
		exit(0);
	}

	conf_params->sol_ID = sol_ID;

	if(sol_ID==SZ)
	{
		layers = (int)iniparser_getint(ini, "PARAMETER:layers", 0);
		conf_params->layers = layers;
		
		int max_quant_intervals = iniparser_getint(ini, "PARAMETER:max_quant_intervals", 65536);
		conf_params->max_quant_intervals = max_quant_intervals;
		maxRangeRadius = max_quant_intervals/2;
		
		stateNum = maxRangeRadius*2;
		allNodes = maxRangeRadius*4;
		
		intvCapacity = maxRangeRadius*2;
		intvRadius = maxRangeRadius;
		
		int quantization_intervals = (int)iniparser_getint(ini, "PARAMETER:quantization_intervals", 0);
		conf_params->quantization_intervals = quantization_intervals;
		if(quantization_intervals>0)
		{
			updateQuantizationInfo(quantization_intervals);
			optQuantMode = 0;
		}
		else
		{
			optQuantMode = 1;
		}
		
		if(quantization_intervals%2!=0)
		{
			printf("Error: quantization_intervals must be an even number!\n");
			iniparser_freedict(ini);
			exit(1);
		}
		
		predThreshold = (float)iniparser_getdouble(ini, "PARAMETER:predThreshold", 0);
		conf_params->predThreshold = predThreshold;
		sampleDistance = (int)iniparser_getint(ini, "PARAMETER:sampleDistance", 0);
		conf_params->sampleDistance = sampleDistance;
		
		offset = (int)iniparser_getint(ini, "PARAMETER:offset", 0);
		conf_params->offset = offset;
		
		modeBuf = iniparser_getstring(ini, "PARAMETER:szMode", NULL);
		if(modeBuf==NULL)
		{
			printf("[SZ] Error: Null szMode setting (please check sz.config file)\n");
			iniparser_freedict(ini);
			exit(1);					
		}
		else if(strcmp(modeBuf, "SZ_BEST_SPEED")==0)
			szMode = SZ_BEST_SPEED;
		else if(strcmp(modeBuf, "SZ_DEFAULT_COMPRESSION")==0)
			szMode = SZ_DEFAULT_COMPRESSION;
		else if(strcmp(modeBuf, "SZ_BEST_COMPRESSION")==0)
			szMode = SZ_BEST_COMPRESSION;
		else
		{
			printf("[SZ] Error: Wrong szMode setting (please check sz.config file)\n");
			iniparser_freedict(ini);
			exit(1);			
		}
		conf_params->szMode = szMode;
		
		modeBuf = iniparser_getstring(ini, "PARAMETER:gzipMode", NULL);
		if(modeBuf==NULL)
		{
			printf("[SZ] Error: Null Gzip mode setting (please check sz.config file)\n");
			iniparser_freedict(ini);
			exit(1);					
		}		
		else if(strcmp(modeBuf, "Gzip_NO_COMPRESSION")==0)
			gzipMode = 0;
		else if(strcmp(modeBuf, "Gzip_BEST_SPEED")==0)
			gzipMode = 1;
		else if(strcmp(modeBuf, "Gzip_BEST_COMPRESSION")==0)
			gzipMode = 9;
		else if(strcmp(modeBuf, "Gzip_DEFAULT_COMPRESSION")==0)
			gzipMode = -1;
		else
		{
			printf("[SZ] Error: Wrong gzip Mode (please check sz.config file)\n");
			exit(0);
		}
		conf_params->gzipMode = gzipMode;
		//maxSegmentNum = (int)iniparser_getint(ini, "PARAMETER:maxSegmentNum", 0); //1024
		
		//spaceFillingCurveTransform = (int)iniparser_getint(ini, "PARAMETER:spaceFillingCurveTransform", 0);
		
		//reOrgSize = (int)iniparser_getint(ini, "PARAMETER:reOrgBlockSize", 0); //8
		
		errBoundMode = iniparser_getstring(ini, "PARAMETER:errorBoundMode", NULL);
		if(errBoundMode==NULL)
		{
			printf("[SZ] Error: Null error bound setting (please check sz.config file)\n");
			iniparser_freedict(ini);
			exit(1);					
		}		
		else if(strcmp(errBoundMode,"ABS")==0||strcmp(errBoundMode,"abs")==0)
			errorBoundMode=ABS;
		else if(strcmp(errBoundMode, "REL")==0||strcmp(errBoundMode,"rel")==0)
			errorBoundMode=REL;
		else if(strcmp(errBoundMode, "ABS_AND_REL")==0||strcmp(errBoundMode, "abs_and_rel")==0)
			errorBoundMode=ABS_AND_REL;
		else if(strcmp(errBoundMode, "ABS_OR_REL")==0||strcmp(errBoundMode, "abs_or_rel")==0)
			errorBoundMode=ABS_OR_REL;
		else
		{
			printf("[SZ] Error: Wrong error bound mode (please check sz.config file)\n");
			iniparser_freedict(ini);
			exit(1);
		}
		conf_params->errorBoundMode = errorBoundMode;
		
		absErrBound = (double)iniparser_getdouble(ini, "PARAMETER:absErrBound", 0);
		conf_params->absErrBound = absErrBound;
		relBoundRatio = (double)iniparser_getdouble(ini, "PARAMETER:relBoundRatio", 0);
		conf_params->relBoundRatio = relBoundRatio;
	}
	
	versionNumber[0] = SZ_VER_MAJOR; //0
	versionNumber[1] = SZ_VER_MINOR; //5
	versionNumber[2] = SZ_VER_REVISION; //15
    
    //initialization for Huffman encoding
    if(pool==NULL)
    {
		pool = (struct node_t*)malloc(allNodes*2*sizeof(struct node_t));
		qqq = (node*)malloc(allNodes*2*sizeof(node));
		code = (unsigned long**)malloc(stateNum*sizeof(unsigned long*));//TODO
		cout = (unsigned char *)malloc(stateNum*sizeof(unsigned char));
		qq = qqq - 1;		
	}
    
    iniparser_freedict(ini);
    return 1;
}

/*-------------------------------------------------------------------------*/
/**
    @brief      It reads and tests the configuration given.
    @return     integer         1 if successfull.

    This function reads the configuration file. Then test that the
    configuration parameters are correct (including directories).

 **/
/*-------------------------------------------------------------------------*/
int SZ_LoadConf() {
    int res = SZ_ReadConf();
    if (res == 0)
    {
        printf("[SZ] ERROR: Impossible to read configuration.\n");
        return 0;
    }
    return 1;
}

int checkVersion(char* version)
{
	int i = 0;
	for(;i<3;i++)
		if(version[i]!=versionNumber[i])
			return 0;
	return 1;
}
