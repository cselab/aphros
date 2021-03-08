#pragma once

#ifdef _MSC_VER
#define VA_SIZE(...) (VA_SIZE0(_0, ## __VA_ARGS__, 2, 1, 0))
#define VA_SIZE0(...) GET_COUNT(__VA_ARGS__)
#define GET_COUNT(_1, _2, _3, N) N
#else
#define GET_COUNT(_1, _2, _3, COUNT, ...) COUNT
#define VA_SIZE(...) GET_COUNT(__VA_ARGS__, 3, 2, 1, 0)
#endif

#define APHROS_CAT(x, y) x##y
#define APHROS_XCAT(x, y) APHROS_CAT(x, y)
#define USEFLAG(x) APHROS_XCAT(0, _USE_##x##_)

#define APHROS_ID_1(x) x
#define APHROS_ID_0(x)

#define MULTIDIMX1 APHROS_XCAT(APHROS_ID_, _USE_DIM1_)(X(1))
#define MULTIDIMX2 APHROS_XCAT(APHROS_ID_, _USE_DIM2_)(X(2))
#define MULTIDIMX3 APHROS_XCAT(APHROS_ID_, _USE_DIM3_)(X(3))
#define MULTIDIMX4 APHROS_XCAT(APHROS_ID_, _USE_DIM4_)(X(4))

#define MULTIDIMX \
  MULTIDIMX1      \
  MULTIDIMX2      \
  MULTIDIMX3      \
  MULTIDIMX4
