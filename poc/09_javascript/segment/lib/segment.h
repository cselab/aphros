#ifdef __cplusplus
#    define EXTERNC extern "C"
#else
#    define EXTERNC
#endif

typedef double Scal;
EXTERNC int segment_get(const Scal*, /**/ Scal **normal, Scal **a);
EXTERNC int segment_norm(int, int, const Scal*, /**/ double*, double*);
