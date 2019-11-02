#define RE 80
#define WE 0.64
// velocity difference over layer
#define SHEAR 0.2
// RHOr = RHO2 / RHO1
#define RHOr 0.01
// MUr = MU2 / MU1
#define MUr 0.01

// liquid
#define RHO1 (1.)
// gas
#define RHO2 (RHO1*RHOr)

#define BCY 0.5
#define BCZ 0.5
#define BR 0.025
#define BRX 1e5
#define BRZ 1e5
#define BCX (1.)
#define BD (BR*2.)

#define SINA 0.1
#define SINF 1
#define SINFZ 1

#define VELC 0.
// characteristic velocity
//#define VEL (VELC)
#define VEL (SHEAR)

// RE = RHO1 * VEL * BD / MU2
// MU2 = RHO1 * VEL * BD / RE
// WE = RHO1 * VEL**2 * BD / SIGMA
// SIGMA = RHO1 * VEL**2 * BD / WE

#define EXTENT 2
#define SQR(a) ((a)*(a))
#define SIGMA (RHO1 * SQR(VEL) * BD / WE)
#define MU1 (RHO1 * VEL * BD / RE)
#define MU2 (MU1*MUr)


#define TMAX (20)
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
