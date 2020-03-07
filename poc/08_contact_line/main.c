#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const char* me = "main";

#define USED(x) \
  if (x)        \
    ;           \
  else {        \
  }

static void usg(void) {
  fprintf(stderr, "echo x y | %s -r ratio -p degree\n", me);
  exit(2);
}

static double sq(double);
static double cube(double);

int main(int argc, char** argv) {
  USED(argc);
  double a;
  double aA;
  double aB;
  double b;
  double bA;
  double bB;
  double c;
  double cA;
  double cB;
  double d;
  double dA;
  double dB;
  double h;
  double mu;
  double pi;
  double psi;
  double r;
  double R;
  double S;
  double C;
  double t;
  double vr;
  double vt;
  double vx;
  double vy;
  double x;
  double y;
  double pressure;
  int Pflag;
  int Rflag;

  Rflag = Pflag = 0;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
      case 'h':
        usg();
        break;
      case 'p':
        argv++;
        if (*argv == NULL) {
          fprintf(stderr, "%s: -p needs an argument\n", me);
          exit(1);
        }
        h = atof(*argv);
        Pflag = 1;
        break;
      case 'r':
        argv++;
        if (*argv == NULL) {
          fprintf(stderr, "%s: -r needs an argument\n", me);
          exit(1);
        }
        R = atof(*argv);
        Rflag = 1;
        break;
      default:
        fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
        exit(2);
    }
  if (!Rflag) {
    fprintf(stderr, "%s: -r is not set\n", me);
    exit(2);
  }
  if (!Pflag) {
    fprintf(stderr, "%s: -p is not set\n", me);
    exit(2);
  }
  if (!(0 <= h && h <= 90)) {
    fprintf(stderr, "%s: h = %g is not in [0, 90]\n", me, h);
    exit(2);
  }

  h *= 0.0174532925199433; /* to radian */
  pi = 3.141592653589793;
  S = sin(h);
  C = cos(h);
  dA = -(S * ((-C * h * pi) + S * pi - C * R * sq(h) + C * sq(h) +
              C * R * sq(S) - C * sq(S))) /
       ((-h * sq(pi)) + C * S * sq(pi) - R * sq(h) * pi + 2 * sq(h) * pi -
        2 * C * S * h * pi + R * sq(S) * pi - C * R * S * sq(h) +
        C * S * sq(h) + R * cube(h) - cube(h) - R * sq(S) * h + sq(S) * h +
        C * R * cube(S) - C * cube(S));
  bA = -dA * pi;
  cA = (C * dA * pi - C * dA * h + S * dA + S) / (S * (h - pi));
  aA = (-cA * pi) - dA - 1;
  bB = 0;
  dB = (dA * h * sq(pi) - 2 * dA * sq(h) * pi + sq(C) * pi - pi + dA * cube(h) +
        sq(C) * dA * h - dA * h) /
       ((sq(h) + sq(C) - 1) * (h - pi));
  aB = (-dB) - 1;
  cB = ((-C * dB * h) + S * dB + S) / (S * h);
  while (scanf("%lf %lf", &x, &y) == 2) {
    t = atan2(y, x);
    r = sqrt(sq(x) + sq(y));
    S = sin(t);
    C = cos(t);
    if (t > h) {
      mu = R;
      a = aA;
      b = bA;
      c = cA;
      d = dA;
    } else {
      mu = 1;
      a = aB;
      b = bB;
      c = cB;
      d = dB;
    }
    vt = (c * t + a) * S + (d * t + b) * C;
    vr = (d * t - c + b) * S - (c * t + d + a) * C;
    psi = r * (c * t * sin(t) + a * sin(t) + d * t * cos(t) + b * cos(t));
    pressure = -(2 * mu * (c * S + d * C)) / r;
    vx = vr * C - vt * S;
    vy = vr * S + vt * C;
    printf("%.16g %.16g %.16g %.16g\n", vx, vy, pressure, psi);
  }
}

static double sq(double x) {
  return x * x;
}

static double cube(double x) {
  return x * x * x;
}
