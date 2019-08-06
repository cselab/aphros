#ifdef CURV_PARTSTR

#ifdef CURV_LEARN
#error Both CURV_PARTSTR and CURV_LEARN are defined
#endif

#include "curvature_partstr.h"

#endif

#ifdef CURV_LEARN
#include "curvature_learn.h"
#endif
