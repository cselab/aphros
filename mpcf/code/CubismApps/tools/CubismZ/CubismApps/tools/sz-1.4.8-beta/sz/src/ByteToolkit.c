/**
 *  @file ByteToolkit.c
 *  @author Sheng Di
 *  @date April, 2016
 *  @brief Byte Toolkit
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */
 
#include <stdlib.h>
#include "sz.h" 
	
int bytesToInt_bigEndian(unsigned char* bytes)
{
	int temp = 0;
	int res = 0;
	
	res <<= 8;
	temp = bytes[0] & 0xff;
	res |= temp;	

	res <<= 8;
	temp = bytes[1] & 0xff;
	res |= temp;

	res <<= 8;
	temp = bytes[2] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = bytes[3] & 0xff;
	res |= temp;		
	
	return res;
}

/**
 * @unsigned char *b the variable to store the converted bytes (length=4)
 * @unsigned int num
 * */
void intToBytes_bigEndian(unsigned char *b, unsigned int num)
{
	b[0] = (unsigned char)(num >> 24);	
	b[1] = (unsigned char)(num >> 16);	
	b[2] = (unsigned char)(num >> 8);	
	b[3] = (unsigned char)(num);	
	
	//note: num >> xxx already considered endian_type...
//if(dataEndianType==LITTLE_ENDIAN_DATA)
//		symTransform_4bytes(*b); //change to BIG_ENDIAN_DATA
}

/**
 * @endianType: refers to the endian_type of unsigned char* b.
 * */
long bytesToLong_bigEndian(unsigned char* b) {
	long temp = 0;
	long res = 0;
	int i;

	res <<= 8;
	temp = b[0] & 0xff;
	res |= temp;

	res <<= 8;
	temp = b[1] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[2] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[3] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[4] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[5] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[6] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[7] & 0xff;
	res |= temp;						
	
	return res;
}

void longToBytes_bigEndian(unsigned char *b, unsigned long num) 
{
	b[0] = (unsigned char)(num>>56);
	b[1] = (unsigned char)(num>>48);
	b[2] = (unsigned char)(num>>40);
	b[3] = (unsigned char)(num>>32);
	b[4] = (unsigned char)(num>>24);
	b[5] = (unsigned char)(num>>16);
	b[6] = (unsigned char)(num>>8);
	b[7] = (unsigned char)(num);
//	if(dataEndianType==LITTLE_ENDIAN_DATA)
//		symTransform_8bytes(*b);
}


long doubleToOSEndianLong(double value)
{
	ldouble buf;
	buf.value = value;
	return buf.lvalue;
}

int floatToOSEndianInt(float value)
{
	lfloat buf;
	buf.value = value;
	return buf.ivalue;
}

//TODO: debug: lfBuf.lvalue could be actually little_endian....
short getExponent_float(float value)
{
	//int ivalue = floatToBigEndianInt(value);

	lfloat lbuf;
	lbuf.value = value;
	int ivalue = lbuf.ivalue;
	
	int expValue = (ivalue & 0x7F800000) >> 23;
	expValue -= 127;
	return (short)expValue;
}

short getPrecisionReqLength_float(float precision)
{
	lfloat lbuf;
	lbuf.value = precision;
	int ivalue = lbuf.ivalue;
	
	int expValue = (ivalue & 0x7F800000) >> 23;
	expValue -= 127;
//	unsigned char the1stManBit = (unsigned char)((ivalue & 0x00400000) >> 22);
//	if(the1stManBit==1)
//		expValue--;	
	return (short)expValue;
}

short getExponent_double(double value)
{
	//long lvalue = doubleToBigEndianLong(value);
	
	ldouble lbuf;
	lbuf.value = value;
	long lvalue = lbuf.lvalue;
	
	int expValue = (int)((lvalue & 0x7FF0000000000000) >> 52);
	expValue -= 1023;
	return (short)expValue;
}

short getPrecisionReqLength_double(double precision)
{
	ldouble lbuf;
	lbuf.value = precision;
	long lvalue = lbuf.lvalue;
	
	int expValue = (int)((lvalue & 0x7FF0000000000000) >> 52);
	expValue -= 1023;
//	unsigned char the1stManBit = (unsigned char)((lvalue & 0x0008000000000000) >> 51);
//	if(the1stManBit==1)
//		expValue--;
	return (short)expValue;
}

unsigned char numberOfLeadingZeros_Int(int i) {
	if (i == 0)
		return 32;
	unsigned char n = 1;
	if (((unsigned int)i) >> 16 == 0) { n += 16; i <<= 16; }
	if (((unsigned int)i) >> 24 == 0) { n +=  8; i <<=  8; }
	if (((unsigned int)i) >> 28 == 0) { n +=  4; i <<=  4; }
	if (((unsigned int)i) >> 30 == 0) { n +=  2; i <<=  2; }
	n -= ((unsigned int)i) >> 31;
	return n;
}

unsigned char numberOfLeadingZeros_Long(long i) {
	 if (i == 0)
		return 64;
	unsigned char n = 1;
	int x = (int)(((unsigned long)i) >> 32);
	if (x == 0) { n += 32; x = (int)i; }
	if (((unsigned int)x) >> 16 == 0) { n += 16; x <<= 16; }
	if (((unsigned int)x) >> 24 == 0) { n +=  8; x <<=  8; }
	if (((unsigned int)x) >> 28 == 0) { n +=  4; x <<=  4; }
	if (((unsigned int)x) >> 30 == 0) { n +=  2; x <<=  2; }
	n -= ((unsigned int)x) >> 31;
	return n;
}

unsigned char getLeadingNumbers_Int(int v1, int v2)
{
	int v = v1 ^ v2;
	return (unsigned char)numberOfLeadingZeros_Int(v);
}

unsigned char getLeadingNumbers_Long(long v1, long v2)
{
	long v = v1 ^ v2;
	return (unsigned char)numberOfLeadingZeros_Long(v);
}

/**
 * By default, the endian type is OS endian type.
 * */
short bytesToShort(unsigned char* bytes)
{
	lshort buf;
	memcpy(buf.byte, bytes, 2);
	
	return buf.svalue;
}

int bytesToInt(unsigned char* bytes)
{
	lfloat buf;
	memcpy(buf.byte, bytes, 4);
	return buf.ivalue;
}

long bytesToLong(unsigned char* bytes)
{
	ldouble buf;
	memcpy(buf.byte, bytes, 8);
	return buf.lvalue;
}

//the byte to input is in the big-endian format
float bytesToFloat(unsigned char* bytes)
{
	lfloat buf;
	memcpy(buf.byte, bytes, 4);
	if(sysEndianType==LITTLE_ENDIAN_SYSTEM)
		symTransform_4bytes(buf.byte);	
	return buf.value;
}

void floatToBytes(unsigned char *b, float num)
{
	lfloat buf;
	buf.value = num;
	memcpy(b, buf.byte, 4);
	if(sysEndianType==LITTLE_ENDIAN_SYSTEM)
		symTransform_4bytes(b);		
}

//the byte to input is in the big-endian format
double bytesToDouble(unsigned char* bytes)
{
	ldouble buf;
	memcpy(buf.byte, bytes, 8);
	if(sysEndianType==LITTLE_ENDIAN_SYSTEM)
		symTransform_8bytes(buf.byte);
	return buf.value;
}

void doubleToBytes(unsigned char *b, double num)
{
	ldouble buf;
	buf.value = num;
	memcpy(b, buf.byte, 8);
	if(sysEndianType==LITTLE_ENDIAN_SYSTEM)
		symTransform_8bytes(b);
}

int extractBytes(unsigned char* byteArray, int k, int validLength)
{
	int outIndex = k/8;
	int innerIndex = k%8;
	unsigned char intBytes[4];
	int length = innerIndex + validLength;
	int byteNum = 0;
	if(length%8==0)
		byteNum = length/8;
	else
		byteNum = length/8+1;
	
	int i;
	for(i = 0;i<byteNum;i++)
		intBytes[4-byteNum+i] = byteArray[outIndex+i];
	
	int result = bytesToInt_bigEndian(intBytes);
	int rightMovSteps = innerIndex +(8 - (innerIndex+validLength)%8)%8;
	result = result << innerIndex;
	switch(byteNum)
	{
	case 1:
		result = result & 0xff;
		break;
	case 2:
		result = result & 0xffff;
		break;
	case 3:
		result = result & 0xffffff;
		break;
	case 4:
		break;
	default: 
		printf("Error: other cases are impossible...\n");
		exit(0);
	}
	result = result >> rightMovSteps;
	
	return result;
}

int getMaskRightCode(int m) {
	switch (m) {
	case 1:
		return 0x01;
	case 2:
		return 0x03;
	case 3:
		return 0x07;
	case 4:
		return 0x0F;
	case 5:
		return 0x1F;
	case 6:
		return 0x3F;
	case 7:
		return 0X7F;
	case 8:
		return 0XFF;
	default:
		return 0;
	}
}

int getLeftMovingCode(int kMod8)
{
	return getMaskRightCode(8 - kMod8);
}

int getRightMovingSteps(int kMod8, int resiBitLength) {
	return 8 - kMod8 - resiBitLength;
}

int getRightMovingCode(int kMod8, int resiBitLength)
{
	int rightMovingSteps = 8 - kMod8 - resiBitLength;
	if(rightMovingSteps < 0)
	{
		switch(-rightMovingSteps)
		{
		case 1:
			return 0x80;
		case 2:
			return 0xC0;
		case 3:
			return 0xE0;
		case 4:
			return 0xF0;
		case 5:
			return 0xF8;
		case 6:
			return 0xFC;
		case 7:
			return 0XFE;
		default:
			return 0;
		}    		
	}
	else //if(rightMovingSteps >= 0)
	{
		int a = getMaskRightCode(8 - kMod8);
		int b = getMaskRightCode(8 - kMod8 - resiBitLength);
		int c = a - b;
		return c;
	}
}

unsigned short* convertByteDataToShortArray(unsigned char* bytes, int byteLength)
{
	lshort ls;
	int i, stateLength = byteLength/2;
	unsigned short* states = (unsigned short*)malloc(stateLength*sizeof(short));
	for(i=0;i<stateLength;i++)
	{
		ls.byte[0] = bytes[i*2];
		ls.byte[1] = bytes[i*2+1];
		states[i] = ls.svalue;
	}
	return states;
} 

void convertShortArrayToBytes(unsigned short* states, int stateLength, unsigned char* bytes)
{
	lshort ls;
	int i, byteLength = byteLength*2;
	for(i=0;i<stateLength;i++)
	{
		ls.svalue = states[i];
		bytes[i*2] = ls.byte[0];
		bytes[i*2+1] = ls.byte[1];
	}
} 
