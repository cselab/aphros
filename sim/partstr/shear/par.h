#define RE 50
#define WE 20
#define RHOr 40
#define MUr 0.5
#define SHEAR 1.

// gas
#define RHO1 (1.)
// liquid
#define RHO2 (RHO1*RHOr)

#define BCY 0.5
#define BCZ 0.5
#define BR 0.2
#define BCX (1.)
#define BD (BR*2.)

#define VELC 0.
#define VEL (BD/SHEAR)

// RE = RHO1 * VEL * BD / MU1
// MU1 = RHO1 * VEL * BD / RE
// WE = RHO1 * VEL**2 * BD / SIGMA
// SIGMA = RHO1 * VEL**2 * BD / WE

#define EXTENT 2
#define SQR(a) ((a)*(a))
#define SIGMA (RHO1 * SQR(VEL) * BD / WE)
#define MU1 (RHO1 * VEL * BD / RE)
#define MU2 (MU1*MUr)

#define TMAX (5)
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
