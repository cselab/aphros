#define REAL double				/* float or double */
#define REALPRINT doubleprint
#define REALRAND doublerand
#define NARROWRAND narrowdoublerand
#define UNIFORMRAND uniformdoublerand

void exactinit(void);
REAL orient2d(REAL *, REAL *, REAL *);
REAL orient2dfast(REAL *, REAL *, REAL *);
REAL orient3d(REAL *, REAL *, REAL *, REAL *);
REAL orient3dfast(REAL *, REAL *, REAL *, REAL *);
REAL incircle(REAL *, REAL *, REAL *, REAL *);
REAL incirclefast(REAL *, REAL *, REAL *, REAL *);
REAL insphere(REAL *, REAL *, REAL *, REAL *, REAL *);
REAL inspherefast(REAL *, REAL *, REAL *, REAL *, REAL *);
