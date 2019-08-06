#ifdef CURV_PARTSTR
  #ifdef CURV_LEARN
  #error Both CURV_PARTSTR and CURV_LEARN are defined
  #endif
  #include "curvature_partstr.h"
#elif defined(CURV_LEARN)
  #include "curvature_learn.h"
#else
  #include "fractions.h"
  #include "curvature.h"
#endif
