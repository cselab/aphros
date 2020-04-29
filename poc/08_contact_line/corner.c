#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const char* me = "corner";

#define USED(x) \
  if (x)        \
    ;           \
  else {        \
  }

static void usg(void) {
  fprintf(stderr, "echo x y | %s -a degree\n", me);
  exit(2);
}

static double sq(double);

int main(int argc, char** argv) {
  USED(argc);
  double a;
  double c;
  double C;
  double D;
  double E;
  double pressure;
  double psi;
  double r;
  double s;
  double S;
  double t;
  double vr;
  double vt;
  double vx;
  double vy;
  double x;
  double X;
  double y;
  int Aflag;
  int mask;

  Aflag = 0;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
      case 'h':
        usg();
        break;
      case 'a':
        argv++;
        if (*argv == NULL) {
          fprintf(stderr, "%s: -a needs an argument\n", me);
          exit(1);
        }
        a = atof(*argv);
        Aflag = 1;
        break;
      default:
        fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
        exit(2);
    }
  if (!Aflag) {
    fprintf(stderr, "%s: -a is not set\n", me);
    exit(2);
  }
  if (!(0 <= a && a <= 180)) {
    fprintf(stderr, "%s: a = %g is not in [0, 180]\n", me, a);
    exit(2);
  }

  a *= 0.0174532925199433; /* to radian */
  C = -sin(a) * sin(a);
  D = -cos(a) * sin(a);
  E = a;
  X = E + D;
  C /= X;
  D /= X;
  E /= X;
  while (scanf("%lf %lf", &x, &y) == 2) {
    t = atan2(y, x);
    mask = t < a;
    if (mask) {
      r = sqrt(sq(x) + sq(y));
      s = sin(t);
      c = cos(t);
      vr = (-D*s*t)+C*c*t+C*s+E*c+D*c;
      vt = (-C*s*t)-D*c*t-E*s;
      psi = r*(C*s*t+D*c*t+E*s);
      pressure = r > 0 ? -2*(C*s+D*c)/r : 0;
      vx = vr * c - vt * s;
      vy = vr * s + vt * c;
      printf("%.16e %.16e %.16e %.16e %d\n", vx, vy, pressure, psi, mask);
    } else
      printf("%.16e %.16e %.16e %.16e %d\n", 0.0, 0.0, 0.0, 0.0, mask);
  }
}

static double sq(double x) {
  return x * x;
}
