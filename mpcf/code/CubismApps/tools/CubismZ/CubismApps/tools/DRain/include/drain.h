#include <stdio.h>

#ifdef __cplusplus
extern "C"
{
#endif
    size_t drain_3df(float * const input,
		     const int xsize,
		     const int ysize,
		     const int zsize,
		     FILE * foutput);

    size_t drain_3df_buffer(float * const input,
		     const int xsize,
		     const int ysize,
		     const int zsize,
		     unsigned char * output);

    size_t pour_3df(FILE * finput,
		    const int xsize,
		    const int ysize,
		    const int zsize,
		    float * const output);

     size_t pour_3df_buffer(unsigned char *input,
		    const int nx,
		    const int ny,
		    const int nz,
		    float * output);

#ifdef __cplusplus
}
#endif
