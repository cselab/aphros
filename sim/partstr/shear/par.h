#define RE 300
#define WE 1000
// RHOr = RHO2 / RHO1
#define RHOr 20
// MUr = MU2 / MU1
#define MUr 20
#define SHEAR 0.05

// gas
#define RHO1 (1.)
// liquid
#define RHO2 (RHO1*RHOr)

#define BCY 0.5
#define BCZ 0.5
#define BR 0.15
#define BCX (BR*1.5)
#define BD (BR*2.)

#define VELC 1.
// characteristic velocity
#define VEL (VELC)
//#define VEL (BR/SHEAR)

// RE = RHO2 * VEL * BD / MU2
// MU2 = RHO2 * VEL * BD / RE
// WE = RHO2 * VEL**2 * BD / SIGMA
// SIGMA = RHO2 * VEL**2 * BD / WE

#define EXTENT 2
#define SQR(a) ((a)*(a))
#define SIGMA (RHO2 * SQR(VEL) * BD / WE)
#define MU2 (RHO2 * VEL * BD / RE)
#define MU1 (MU2/MUr)

#define TMAX (10)
#define DUMPDT (0.1)

#ifndef NX
#define NX 64
#endif


#ifndef MPIDIM
#define MPIDIM 
#endif

#ifndef NOIO
#define NOIO 0
#endif
