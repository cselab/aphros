#pragma once

#define GET_COUNT(_1, _2, _3, COUNT, ...) COUNT
#define VA_SIZE(...) GET_COUNT(__VA_ARGS__, 3, 2, 1, 0)
#define APHROS_CAT(x, y) x##y
#define APHROS_XCAT(x, y) APHROS_CAT(x, y)
#define USEFLAG(x) APHROS_XCAT(0, _USE_##x##_)

