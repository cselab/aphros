#pragma once

#ifdef __cplusplus
extern "C"
{
#endif
    void iw4_gen2_fw(float * data,
			  const int xsize,
			  const int ysize,
			  const int zsize,
			  int& xsizecoarse,
			  int& ysizecoarse,
			  int& zsizecoarse);
    
    void iw4_gen2_bw(float * data,
			 const int xsize,
			 const int ysize,
			 const int zsize,
			 const int xsizecoarse,
			 const int ysizecoarse,
			 const int zsizecoarse);
#ifdef __cplusplus
}
#endif
