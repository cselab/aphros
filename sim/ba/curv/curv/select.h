#ifdef CURV_PARTSTR
  #warning Using partstr.h
  #include "partstr.h"
#endif

#if LEARN_N10
  #warning model sw3d3_r32_s322_n10
  #include "learn_eval_sw3d3_r32_s322_n10.h"
#endif

#if LEARN_N80
  #warning model sw3d3_r32_s322_n80
  #include "learn_eval_sw3d3_r32_s322_n80.h"
#endif

#if LEARN_N80_SPHCYL
  #warning model sw3d3sphcyl_r32_s440
  #include "learn_eval_sw3d3sphcyl_r32_s440.h"
#endif

#ifdef CURV_LEARN
  #warning Using learn.h
  #include "learn.h"
#endif

#ifdef CURV_DIV
  #warning Using div.h
  #include "div.h"
#endif

#ifdef CURV_FIX
  #warning Using fix.h
  #include "fix.h"
#endif

#include "fractions.h"
#include "curvature.h"
