#define RE 50
#define WE 50
#define RHOr 20
#define MUr 1
#define SHEAR 0.625

// gas
#define RHO1 (1.)
// liquid
#define RHO2 (RHO1*RHOr)

#define BCY 0.5
#define BCZ 0.5
#define BR 0.15
#define BCX (BR*2)
#define BD (BR*2)

#define VEL 1.

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
