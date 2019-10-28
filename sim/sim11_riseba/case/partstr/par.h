#define NX 256

#define GA 1000
#define EO 20

#define RHO1 (1.)
#define RHO2 (RHO1*0.01)

#define BCY 0.5
#define BCZ 0.5
#define BR 0.15
#define BCX (BR*2)

// Ga = (rho0 * sqrt(g*R) * R / mu0)
// g = (Ga * mu0 / (rho0 * R))**2 / R
// mu0 = rho0 * sqrt(g*R) * R  / Ga
// Eo = (rho0 * g * R*R / sigma)
// sigma = (rho0 * g * R*R) / (Eo)

#define EXTENT 2
#define SQR(a) ((a)*(a))
#define GRAV 1
#define SIGMA ((RHO1 * GRAV * BR * BR) / (EO))
#define MU1 (RHO1 * sqrt(GRAV * BR) * BR / GA)
#define MU2 (MU1*0.01)

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
