#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "myzopfli.h"
#include <zlib.h>

#include <omp.h>

#ifdef _USE_DEFLATE_
int CompressData( const char* abSrc, int nLenSrc, char* abDst, int nLenDst )
{
    z_stream zInfo ={0};
    zInfo.total_in=  zInfo.avail_in=  nLenSrc;
    zInfo.total_out= zInfo.avail_out= nLenDst;
    zInfo.next_in= (char*)abSrc;
    zInfo.next_out= abDst;

    int nErr, nRet= -1;
//    nErr= deflateInit( &zInfo, Z_DEFAULT_COMPRESSION ); // zlib function
    nErr= deflateInit( &zInfo, Z_BEST_COMPRESSION ); // zlib function
    if ( nErr == Z_OK ) {
        nErr= deflate( &zInfo, Z_FINISH );              // zlib function
        if ( nErr == Z_STREAM_END ) {
            nRet= zInfo.total_out;
        }
    }
    deflateEnd( &zInfo );    // zlib function
    return( nRet );
}

int UncompressData( const char* abSrc, int nLenSrc, char* abDst, int nLenDst )
{
    z_stream zInfo ={0};
    zInfo.total_in=  zInfo.avail_in=  nLenSrc;
    zInfo.total_out= zInfo.avail_out= nLenDst;
    zInfo.next_in= (char*)abSrc;
    zInfo.next_out= abDst;

    int nErr, nRet= -1;
    nErr= inflateInit( &zInfo );               // zlib function
    if ( nErr == Z_OK ) {
        nErr= inflate( &zInfo, Z_FINISH );     // zlib function
        if ( nErr == Z_STREAM_END ) {
            nRet= zInfo.total_out;
        }
    }
    inflateEnd( &zInfo );   // zlib function
    return( nRet ); // -1 or len of output
}
#endif 

void shuffle(char *in, int n, int s)
{
	int i, j, k;
	int b;

	char *tmp = malloc(n);

	j = 0;
	for (b = 0; b < s; b++)
	{
		for (i = b; i < n; i+= s) {
			tmp[j] = in[i];
			j++;
		}
	}

	memcpy(in, tmp, n);
	free(tmp);
}

void reshuffle(char *in, int n, int s)
{
	int i, j, k;
	int b;

	char *tmp = malloc(n);
	int c = n/s;
	
	j = 0;
	for (i = 0; i < n/s; i++) {
		for (b = 0; b < s; b++) {
			tmp[j] = in[i + b*c]; j++;

		}
	}

	memcpy(in, tmp, n);
	free(tmp);
}


#ifndef dtype
#define dtype float
#endif

#define SZ	(32*32*32)

static dtype indat[SZ];
static dtype outdat[SZ];
static dtype outdat2[SZ];

int main(int argc, char *argv[])
{
//	unsigned char *myout = NULL;
        size_t outs;
        size_t outs2;

        int i;
//	for (i = 0; i < SZ; i++) indat[i] = i;
	srand48(1);
	for (i = 0; i < SZ; i++) indat[i] = drand48();

//	shuffle((char *)indat, SZ*sizeof(dtype), sizeof(dtype));


	double t0 = omp_get_wtime();
#ifndef _USE_DEFLATE_
	outs = zopfli_compress((const unsigned char *) indat, SZ*sizeof(dtype), (unsigned char *) outdat);
#else
	outs = CompressData( (const char*)indat, SZ*sizeof(dtype), (char *) outdat, 2*SZ*sizeof(dtype));
#endif
	double t1 = omp_get_wtime();
	
        printf("outs = %ld\n", outs);
        printf("rate = %.2lf\n", (1.0*SZ*sizeof(dtype))/outs);
        printf("time = %.2lf s\n", t1-t0);

#ifndef _USE_DEFLATE_
	outs2 = zopfli_decompress((const void *)outdat, outs, (void *) outdat2, SZ*sizeof(dtype));
#else
	outs2 = zopfli_decompress((const void *)outdat, outs, (void *) outdat2, SZ*sizeof(dtype));
//	outs2 = UncompressData( (const char*) outdat, outs, (char *) outdat2, SZ*sizeof(dtype));
#endif
//	outs2 = SZ*sizeof(dtype);

        printf("outs2 = %d\n", outs2);
//	reshuffle((char *)outdat2, SZ*sizeof(dtype), sizeof(dtype));

        for (i = 0; i < SZ; i+=(SZ/4)) printf("%d : %lf\n", i, outdat2[i]);
        for (i = SZ-3; i < SZ; i+=1) printf("%d : %lf\n", i, outdat2[i]);

 return 0;
}

