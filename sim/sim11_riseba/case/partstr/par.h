#define NX 128

#define EO 20
#define GA 70.71

#define RHO1 (1000.)
#define RHO2 (1.)
#define MU1 (10.)
#define MU2 (0.1)

#define BCX 0.3
#define BCY 0.3
#define BCZ 0.3
#define BR 0.2

// Ga = (rho0 * sqrt(g*R) * R / mu0)
// g = (Ga * mu0 / (rho0 * R))**2 / R
// Eo = (rho0 * g * R*R / sigma)
// sigma = (rho0 * g * R*R) / (Eo)

#define TMAX (3.)
#define DUMPDT (0.1)
#define EXTENT 2.
#define SQR(a) ((a)*(a))
#define GRAV (SQR(GA * MU1 / (RHO1 * BR)) / BR)
#define SIGMA ((RHO1 * GRAV * BR * BR) / (EO))

#ifndef NX
#define NX 64
#endif

//#define MPIDIM dimensions(nx = 4, ny = 4, nz = 4);

#ifndef MPIDIM
#define MPIDIM 
#endif

#ifndef NOIO
#define NOIO 0
#endif
