#define TMAX (10.0)
#define DUMPDT (0.5)
#define SIGMA (0.1)
#define RHO1 (1.0)
#define RHO2 (1.0)
#define MU1 (0.000625)
#define MU2 (0.000625)
#define EXTENT (6.283185307179586)

#define BCX 2.
#define BCY 2.
#define BCZ 2.
#define BR 0.2

#ifndef NX
#define NX 64
#endif

//#define MPIDIM dimensions(nx = 4, ny = 4, nz = 4);

#ifndef MPIDIM
#define MPIDIM 
#endif
