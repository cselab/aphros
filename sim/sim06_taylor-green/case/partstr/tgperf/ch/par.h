// generated by "./runcase nx=256 np=512 partstr=1 suff=_fixdt"
#define NX 256
#define NOIO 0
#define CURV_PARTSTR
#define TMAX (10.0)
#define DUMPDT (1.0)
#define SIGMA (0.1)
#define RHO1 (1.0)
#define RHO2 (0.01)
#define MU1 (0.00125)
#define MU2 (0.0000125)
#define EXTENT (6.283185307179586)

#define BCX 2.
#define BCY 2.
#define BCZ 2.
#define BR 0.1

#ifndef NX
#define NX 64
#endif

//#define MPIDIM dimensions(nx = 4, ny = 4, nz = 4);

#define MPIDIM DT=0.001;
//#define MPIDIM DT=0.003; kPartstr.Np = 5;  kPartstr.Ns = 2; kPartstr.eps = 1e-4;

#ifndef MPIDIM
#define MPIDIM 
#endif

#ifndef NOIO
#define NOIO 0
#endif

#define BDAT "../b890.dat"
