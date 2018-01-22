/* 
Copyright © 2006 Cornell Research Foundation, Inc.  All rights reserved.
Author: Professor Martin Burtscher

Software License Terms and Conditions
1. SOFTWARE shall mean the FPC code, or portions thereof, available via the web page http://www.csl.cornell.edu/~burtscher/research/FPC/ and described in Cornell Research Foundation, Inc. („CRF‰) file D-3944.  SOFTWARE includes, but is not limited to, source code, object code and executable code.
2. CRF is a wholly owned subsidiary of Cornell University, is a fiduciary of Cornell University in intellectual property matters and holds all intellectual property rights in SOFTWARE.
3. LICENSEE means the party to this Agreement and the user of SOFTWARE.  By using SOFTWARE, LICENSEE enters into this Agreement with CRF.
4. SOFTWARE is made available under this Agreement to allow certain non-commercial research and teaching use.  CRF reserves all commercial rights to SOFTWARE and these rights may be licensed by CRF to third parties.
5. LICENSEE is hereby granted permission to: a) use SOFTWARE for non-commercial research or teaching purposes, and b) download, compile, execute, copy, and modify SOFTWARE for non-commercial research or teaching purposes provided that this notice accompanies all copies of SOFTWARE.  Copies of modified SOFTWARE may be distributed only for non-commercial research or teaching purposes (i) if this notice accompanies those copies, (ii) if said copies carry prominent notices stating that SOFTWARE has been changed, and (iii) the date of any changes are clearly identified in SOFTWARE.
6. CRF may terminate this Agreement at any time if LICENSEE breaches a material provision of this Agreement.  CRF may also terminate this Agreement if the SOFTWARE becomes subject to any claim of infringement of patent, copyright or trade secret, or if in CRF‚s opinion such a claim is likely to occur.
7. LICENSEE agrees that the export of SOFTWARE from the United States may require approval from the U.S. government and failure to obtain such approval will result in the immediate termination of this license and may result in criminal liability under U.S. laws.
8. The work leading to the development of SOFTWARE was supported in part by various grants from an agency of the U.S. Government, and CRF is obligated to comply with U.S. OMB Circular A-124 and 37 CFR Part 401.  This license is subject to the applicable terms of U.S. Government regulations concerning Government funded inventions.
9. CRF provides SOFTWARE on an „as is‰ basis.  CRF does not warrant, guarantee, or make any representations regarding the use or results of SOFTWARE with respect to its correctness, accuracy, reliability or performance.  The entire risk of the use and performance of SOFTWARE is assumed by LICENSEE.  ALL WARRANTIES INCLUDING, WITHOUT LIMITATION, ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE OR MERCHANTABILITY AND ANY WARRANTY OF NONINFRINGEMENT OF PATENTS, COPYRIGHTS, OR ANY OTHER INTELLECTUAL PROPERTY RIGHT ARE HEREBY EXCLUDED.
10. LICENSEE understands and agrees that neither CRF nor Cornell University is under any obligation to provide maintenance, support or update services, notices of latent defects, correction of defects, or future versions for SOFTWARE.
11. Even if advised of the possibility of damages, under no circumstances shall CRF or Cornell University individually or jointly be liable to LICENSEE or any third party for damages of any character, including, without limitation, direct, indirect, incidental, consequential or special damages, loss of profits, loss of use, loss of goodwill, computer failure or malfunction.  LICENSEE agrees to indemnify and hold harmless CRF and Cornell University for any and all liability CRF or Cornell University may incur as a result of use of SOFTWARE by LICENSEE.
*/


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#define SIZE 32768

static const long long mask[8] =
{0x0000000000000000LL,
0x00000000000000ffLL,
0x000000000000ffffLL,
0x0000000000ffffffLL,
0x000000ffffffffffLL,
0x0000ffffffffffffLL,
0x00ffffffffffffffLL,
0xffffffffffffffffLL};

void fpc_compress(char datain[], long inbytes, char dataout[], int *outbytes, long predsizem1)
{
 register long i, out, intot, hash, dhash, code, bcode, ioc;
 register long long val, lastval, stride, pred1, pred2, xor1, xor2;
 register long long *fcm, *dfcm;
 unsigned long long inbuf[SIZE + 1];
 unsigned char outbuf[6 + (SIZE / 2) + (SIZE * 8) + 2];

 int bytes_written = 0;
 int bytes_read = 0;
 int bytes_rem;

 assert(0 == ((long)outbuf & 0x7));

 outbuf[0] = predsizem1;
#if 0
 ioc = fwrite(outbuf, 1, 1, stdout);
 assert(1 == ioc);
#else
 memcpy(&dataout[bytes_written], outbuf, 1); bytes_written+=1;
#endif
 predsizem1 = (1L << predsizem1) - 1;

 hash = 0;
 dhash = 0;
 lastval = 0;
 pred1 = 0;
 pred2 = 0;
 fcm = (long long *)calloc(predsizem1 + 1, 8);
 assert(NULL != fcm);
 dfcm = (long long *)calloc(predsizem1 + 1, 8);
 assert(NULL != dfcm);

#if 0
 intot = fread(inbuf, 8, SIZE, stdin);
#else
 bytes_rem = inbytes - bytes_read;
 if (bytes_rem > 0) {
	intot = (bytes_rem < 8*SIZE) ? bytes_rem : 8*SIZE;
	memcpy(inbuf, datain+bytes_read, intot);
	bytes_read += intot;
	intot = intot/8;
 }
 else {
	 intot = 0;
 }
#endif
 while (0 < intot) {
   val = inbuf[0];
   out = 6 + ((intot + 1) >> 1);
   *((long long *)&outbuf[(out >> 3) << 3]) = 0;
   for (i = 0; i < intot; i += 2) {
     xor1 = val ^ pred1;
     fcm[hash] = val;
     hash = ((hash << 6) ^ ((unsigned long long)val >> 48)) & predsizem1;
     pred1 = fcm[hash];

     stride = val - lastval;
     xor2 = val ^ (lastval + pred2);
     lastval = val;
     val = inbuf[i + 1];
     dfcm[dhash] = stride;
     dhash = ((dhash << 2) ^ ((unsigned long long)stride >> 40)) & predsizem1;
     pred2 = dfcm[dhash];

     code = 0;
     if ((unsigned long long)xor1 > (unsigned long long)xor2) {
       code = 0x80;
       xor1 = xor2;
     }
     bcode = 7;                // 8 bytes
     if (0 == (xor1 >> 56))
       bcode = 6;              // 7 bytes
     if (0 == (xor1 >> 48))
       bcode = 5;              // 6 bytes
     if (0 == (xor1 >> 40))
       bcode = 4;              // 5 bytes
     if (0 == (xor1 >> 24))
       bcode = 3;              // 3 bytes
     if (0 == (xor1 >> 16))
       bcode = 2;              // 2 bytes
     if (0 == (xor1 >> 8))
       bcode = 1;              // 1 byte
     if (0 == xor1)
       bcode = 0;              // 0 bytes

     *((long long *)&outbuf[(out >> 3) << 3]) |= xor1 << ((out & 0x7) << 3);
     if (0 == (out & 0x7))
       xor1 = 0;
     *((long long *)&outbuf[((out >> 3) << 3) + 8]) = (unsigned long long)xor1 >> (64 - ((out & 0x7) << 3));

     out += bcode + (bcode >> 2);
     code |= bcode << 4;

     xor1 = val ^ pred1;
     fcm[hash] = val;
     hash = ((hash << 6) ^ ((unsigned long long)val >> 48)) & predsizem1;
     pred1 = fcm[hash];

     stride = val - lastval;
     xor2 = val ^ (lastval + pred2);
     lastval = val;
     val = inbuf[i + 2];
     dfcm[dhash] = stride;
     dhash = ((dhash << 2) ^ ((unsigned long long)stride >> 40)) & predsizem1;
     pred2 = dfcm[dhash];

     bcode = code | 0x8;
     if ((unsigned long long)xor1 > (unsigned long long)xor2) {
       code = bcode;
       xor1 = xor2;
     }
     bcode = 7;                // 8 bytes
     if (0 == (xor1 >> 56))
       bcode = 6;              // 7 bytes
     if (0 == (xor1 >> 48))
       bcode = 5;              // 6 bytes
     if (0 == (xor1 >> 40))
       bcode = 4;              // 5 bytes
     if (0 == (xor1 >> 24))
       bcode = 3;              // 3 bytes
     if (0 == (xor1 >> 16))
       bcode = 2;              // 2 bytes
     if (0 == (xor1 >> 8))
       bcode = 1;              // 1 byte
     if (0 == xor1)
       bcode = 0;              // 0 bytes

     *((long long *)&outbuf[(out >> 3) << 3]) |= xor1 << ((out & 0x7) << 3);
     if (0 == (out & 0x7))
       xor1 = 0;
     *((long long *)&outbuf[((out >> 3) << 3) + 8]) = (unsigned long long)xor1 >> (64 - ((out & 0x7) << 3));

     out += bcode + (bcode >> 2);
     outbuf[6 + (i >> 1)] = code | bcode;
   }
   if (0 != (intot & 1)) {
     out -= bcode + (bcode >> 2);
   }
   outbuf[0] = intot;
   outbuf[1] = intot >> 8;
   outbuf[2] = intot >> 16;
   outbuf[3] = out;
   outbuf[4] = out >> 8;
   outbuf[5] = out >> 16;
#if 0
   ioc = fwrite(outbuf, 1, out, stdout);
   assert(ioc == out);
#else
	memcpy(&dataout[bytes_written], outbuf, out); bytes_written+=out;
#endif

#if 0
   intot = fread(inbuf, 8, SIZE, stdin);
#else
	bytes_rem = inbytes - bytes_read;
	if (bytes_rem > 0) {
		intot = (bytes_rem < 8*SIZE) ? bytes_rem : 8*SIZE;
		memcpy(inbuf, datain+bytes_read, intot);
		bytes_read += intot;
		intot = intot/8;
	}
	else {
		intot = 0;
	}
#endif
 }
 
 *outbytes = bytes_written;
 
 return;
}

void fpc_decompress(char datain[], long inbytes, char dataout[], int *outbytes)
{
 register long i, in, intot, hash, dhash, code, bcode, predsizem1, end, tmp, ioc;
 register long long val, lastval, stride, pred1, pred2, next;
 register long long *fcm, *dfcm;
 long long outbuf[SIZE];
 unsigned char inbuf[(SIZE / 2) + (SIZE * 8) + 6 + 2];

 int bytes_written = 0;
 int bytes_read = 0;
 int bytes_rem;

 assert(0 == ((long)inbuf & 0x7));

#if 0
 ioc = fread(inbuf, 1, 7, stdin);
#else
 bytes_rem = inbytes - bytes_read;
 if (bytes_rem > 0) {
	memcpy(inbuf, datain+bytes_read, 7);
	bytes_read += 7;
	ioc = 7;
 }
#endif
 if (1 != ioc) {
   assert(7 == ioc);
   predsizem1 = inbuf[0];
   predsizem1 = (1L << predsizem1) - 1;

   hash = 0;
   dhash = 0;
   lastval = 0;
   pred1 = 0;
   pred2 = 0;
   fcm = (long long *)calloc(predsizem1 + 1, 8);
   assert(NULL != fcm);
   dfcm = (long long *)calloc(predsizem1 + 1, 8);
   assert(NULL != dfcm);

   intot = inbuf[3];
   intot = (intot << 8) | inbuf[2];
   intot = (intot << 8) | inbuf[1];
   in = inbuf[6];
   in = (in << 8) | inbuf[5];
   in = (in << 8) | inbuf[4];
   assert(SIZE >= intot);
   do {
#if 0
		end = fread(inbuf, 1, in, stdin);
#else
		bytes_rem = inbytes - bytes_read;
		if (bytes_rem > 0) {
			end = (bytes_rem < in) ? bytes_rem : in;
			memcpy(inbuf, datain+bytes_read, end);
			bytes_read += end;
		}
		else {
			printf("critical error!\n");
			exit(1);
		}
#endif
     assert((end + 6) >= in);
     in = (intot + 1) >> 1;
     for (i = 0; i < intot; i += 2) {
       code = inbuf[i >> 1];

       val = *((long long *)&inbuf[(in >> 3) << 3]);
       next = *((long long *)&inbuf[((in >> 3) << 3) + 8]);
       tmp = (in & 0x7) << 3;
       val = (unsigned long long)val >> tmp;
       next <<= 64 - tmp;
       if (0 == tmp)
         next = 0;
       val |= next;

       bcode = (code >> 4) & 0x7;
       val &= mask[bcode];
       in += bcode + (bcode >> 2);

       if (0 != (code & 0x80))
         pred1 = pred2;
       val ^= pred1;

       fcm[hash] = val;
       hash = ((hash << 6) ^ ((unsigned long long)val >> 48)) & predsizem1;
       pred1 = fcm[hash];

       stride = val - lastval;
       dfcm[dhash] = stride;
       dhash = ((dhash << 2) ^ ((unsigned long long)stride >> 40)) & predsizem1;
       pred2 = val + dfcm[dhash];
       lastval = val;

       outbuf[i] = val;

       val = *((long long *)&inbuf[(in >> 3) << 3]);
       next = *((long long *)&inbuf[((in >> 3) << 3) + 8]);
       tmp = (in & 0x7) << 3;
       val = (unsigned long long)val >> tmp;
       next <<= 64 - tmp;
       if (0 == tmp)
         next = 0;
       val |= next;

       bcode = code & 0x7;
       val &= mask[bcode];
       in += bcode + (bcode >> 2);

       if (0 != (code & 0x8))
         pred1 = pred2;
       val ^= pred1;

       fcm[hash] = val;
       hash = ((hash << 6) ^ ((unsigned long long)val >> 48)) & predsizem1;
       pred1 = fcm[hash];

       stride = val - lastval;
       dfcm[dhash] = stride;
       dhash = ((dhash << 2) ^ ((unsigned long long)stride >> 40)) & predsizem1;
       pred2 = val + dfcm[dhash];
       lastval = val;

       outbuf[i + 1] = val;
     }
#if 0
     ioc = fwrite(outbuf, 8, intot, stdout);
     assert(ioc == intot);
#else
	memcpy(&dataout[bytes_written], outbuf, 8*intot); bytes_written+=8*intot;
#endif
     intot = 0;
     if ((end - 6) >= in) {
       intot = inbuf[in + 2];
       intot = (intot << 8) | inbuf[in + 1];
       intot = (intot << 8) | inbuf[in];
       end = inbuf[in + 5];
       end = (end << 8) | inbuf[in + 4];
       end = (end << 8) | inbuf[in + 3];
       in = end;
     }
     assert(SIZE >= intot);
   } while (0 < intot);
 }
 *outbytes = bytes_written;
}

