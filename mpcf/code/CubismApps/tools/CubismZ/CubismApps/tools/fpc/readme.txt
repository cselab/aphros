FPC is a very fast and effective lossless compressor/decompressor for IEEE 754 64-bit double-precision floating-point data.

Enter your email address below (which will be stored in our database) and press the submit button to have the source code emailed to you.
A description of FPC is available here.  A line-by-line description of the source code is available here.
Sample little-endian datasets are available here.  Note that FPC is protected by this license.

The source code in file fpc.c can be compiled into an executable called fpc as follows:

gcc -O3 fpc.c -o fpc

The executable compresses/decompresses data with level x from the standard input and writes the result to the standard output.
x is a positive integer that specifies the internal table size (2x+4 bytes).
Larger x typically result in better compression ratios but slower processing.

To compress the file file.in with level 20 and store the result in a file called file.fpc, enter:

./fpc 20 < file.in > file.fpc

To decompress the file file.fpc and store the result in file file.out, enter:

./fpc < file.fpc > file.out

Note that the raw data file has to be a multiple of 8 bytes long and should contain nothing but binary double-precision values.
Only little-endian systems are currently supported.
FPC works particularly well on 64-bit machines and on CPUs with at least 16 general-purpose integer registers.
