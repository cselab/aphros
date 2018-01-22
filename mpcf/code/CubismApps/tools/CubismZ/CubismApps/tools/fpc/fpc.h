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


void fpc_compress(char datain[], long inbytes, char dataout[], int *outbytes, long predsizem1);
void fpc_decompress(char datain[], long inbytes, char dataout[], int *outbytes);
