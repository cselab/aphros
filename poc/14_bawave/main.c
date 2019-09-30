#include ".u/io/io.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"

#define DUMPDT (0.01)
#define TMAX (1)
#define ak  (0.55)
#define h_   (0.5)
#define k_  (2.*pi)

static double wave(double, double);
int main() {
  init_grid(64);
  origin (-L0/2, -L0/2, -L0/2);
  periodic (right);

  rho1 = 1;
  rho2 = 0.01;

  mu1 = 0;
  mu2 = 0;

  run();
}

event init (i = 0) {
  foreach() {
    u.x[] = 1;
    u.y[] = 0;
  }
  fraction(f, wave(x, y));
}

event out (t += DUMPDT ; t <= TMAX + DUMPDT) {
  fprintf(stderr, "dump i=%05d t=%g dt=%g \n", i, t ,dt);

  static int frame = 0;
  scalar * a = {u, p, f};

  char name[1000];

  sprintf(name, "u_%04d.vtk", frame);
  FILE * fp = fopen(name, "w");
  io(a, fp);
  fclose(fp);

  ++frame;
}

event statout (i += 10 ; t <= TMAX) {
  fprintf(stderr, "step i=%05d t=%g dt=%g \n", i, t ,dt);
}

static double eta (double x, double y)
{
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + ak*eta2 + sq(ak)*eta3;
}

static double
wave0(double x, double y)
{
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + ak*eta2 + sq(ak)*eta3 - y;
}

static double
wave(double x, double y)
{
    return wave0(x, y);
}
