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
  double psi;
  double x;
  double y;
  double vx;
  double vy;
  double vt;
  double vr;
  double D;
  double E;
  double S;
  double C;
  double a;
  double t;
  double r;
  double X;
  double pressure;
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
  while (scanf("%lf %lf", &x, &y) == 2) {
    t = atan2(y, x);
    mask = t < a;
    if (mask) {
      r = sqrt(sq(x) + sq(y));
      S = sin(t);
      C = cos(t);
      vr = ((-D * S * t) + C * C * t + C * S + C * E + C * D) / X;
      vt = -(C * S * t + C * D * t + E * S) / X;
      psi = (r * (C * S * t + C * D * t + E * S)) / X;
      pressure = r > 0 ? -(2 * (C * sin(t) + D * cos(t))) / (X * r) : 0;
      vx = vr * C - vt * S;
      vy = vr * S + vt * C;
      printf("%.16e %.16e %.16e %.16e %d\n", vx, vy, pressure, psi, mask);
    } else
      printf("%.16e %.16e %.16e %.16e %d\n", 0.0, 0.0, 0.0, 0.0, mask);
  }
}

static double sq(double x) {
  return x * x;
}
