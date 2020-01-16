// Created by Sergey Litvinov on 28.09.2019
// Copyright 2019 ETH Zurich

#include <stdio.h>
#include <tgmath.h>

static double sq(double x) {
  return x * x;
}

#define pi (3.141592653589793)

static double k = 2 * pi;
static double h = 0.5;
static double g = 1;
static double a = 0.1;
static double ak = 0.6283185307179586;

void f0(double x, double y, double* pvx, double* pvy, double* peta) {
  double e, c, o;
  double vx, vy, eta;

  e = a * k;
  c = 1 / tanh(h * k);
  eta =
      a * sq(e) *
          ((3.0E+0 * (sq(sq(c) - 1) + 8 * pow(c, 6)) * cos(3 * k * x)) /
               6.4E+1 +
           ((-3.0E+0) * ((-3 * sq(c)) + pow(c, 4) + 3) * cos(k * x)) / 8.0E+0) +
      (a * c * (3 * sq(c) - 1) * e * cos(2 * k * x)) / 4.0E+0 + a * cos(k * x);
  o = sqrt(g * k * tanh(h * k)) *
      sqrt(
          sq(a) * sq(k) *
              ((9.0E+0 * sq(1 / pow(tanh(h * k), 2) - 1)) / 8.0E+0 +
               1 / pow(tanh(h * k), 2)) +
          1);
  vx = (3.0E+0 * a * (sq(c) - 1) * (sq(c) + 3) * (9 * sq(c) - 13) * sq(e) * g *
        k * cos(3 * k * x) * cosh(3 * k * (y + h))) /
           (6.4E+1 * cosh(3 * h * k) * o) +
       (3.0E+0 * a * sq(sq(c) - 1) * e * g * k * cos(2 * k * x) *
        cosh(2 * k * (y + h))) /
           (4.0E+0 * c * cosh(2 * h * k) * o) +
       (a * g * k * cos(k * x) * cosh(k * (y + h))) / (cosh(h * k) * o);
  vy = (3.0E+0 * a * (sq(c) - 1) * (sq(c) + 3) * (9 * sq(c) - 13) * sq(e) * g *
        k * sin(3 * k * x) * sinh(3 * k * (y + h))) /
           (6.4E+1 * cosh(3 * h * k) * o) +
       (3.0E+0 * a * sq(sq(c) - 1) * e * g * k * sin(2 * k * x) *
        sinh(2 * k * (y + h))) /
           (4.0E+0 * c * cosh(2 * h * k) * o) +
       (a * g * k * sin(k * x) * sinh(k * (y + h))) / (cosh(h * k) * o);

  *pvx = vx;
  *pvy = vy;
  *peta = eta;
}

void f1(double x, double y, double* pvx, double* pvy, double* peta) {
  double eps, chi, omega;
  double vx, vy, eta;

  eps = a * k;

  chi = 1.0 / tanh(h * k);

  eta = (1.0L / 4.0L) * a * chi * eps * (3 * pow(chi, 2) - 1) * cos(2 * k * x) +
        a * pow(eps, 2) *
            ((1.0L / 64.0L) * (24 * pow(chi, 6) + 3 * pow(pow(chi, 2) - 1, 2)) *
                 cos(3 * k * x) +
             (1.0L / 8.0L) * (-3 * pow(chi, 4) + 9 * pow(chi, 2) - 9) *
                 cos(k * x)) +
        a * cos(k * x);

  omega = sqrt(
      g * k *
      (pow(eps, 2) * (pow(chi, 2) + (9.0L / 8.0L) * pow(pow(chi, 2) - 1, 2)) +
       1) *
      tanh(h * k));

  vx = (3.0L / 64.0L) * a * pow(eps, 2) * g * k * (pow(chi, 2) - 1) *
           (pow(chi, 2) + 3) * (9 * pow(chi, 2) - 13) * cos(3 * k * x) *
           cosh(3 * k * (h + y)) / (omega * cosh(3 * h * k)) +
       a * g * k * cos(k * x) * cosh(k * (h + y)) / (omega * cosh(h * k)) +
       (3.0L / 4.0L) * a * eps * g * k * pow(pow(chi, 2) - 1, 2) *
           cos(2 * k * x) * cosh(2 * k * (h + y)) /
           (chi * omega * cosh(2 * h * k));

  vy = (3.0L / 64.0L) * a * pow(eps, 2) * g * k * (pow(chi, 2) - 1) *
           (pow(chi, 2) + 3) * (9 * pow(chi, 2) - 13) * sin(3 * k * x) *
           sinh(3 * k * (h + y)) / (omega * cosh(3 * h * k)) +
       a * g * k * sin(k * x) * sinh(k * (h + y)) / (omega * cosh(h * k)) +
       (3.0L / 4.0L) * a * eps * g * k * pow(pow(chi, 2) - 1, 2) *
           sin(2 * k * x) * sinh(2 * k * (h + y)) /
           (chi * omega * cosh(2 * h * k));

  *pvx = vx;
  *pvy = vy;
  *peta = eta;
}

void f2(double x, double y, double* p) {
  double Phi;
  double alpa = 1. / tanh(k * h);
  double a = ak / k;
  double sgma = sqrt(
      g * k * tanh(k * h) *
      (1. + k * k * a * a *
                (9. / 8. * (sq(alpa) - 1.) * (sq(alpa) - 1.) + sq(alpa))));
  double A = a * g / sgma;
  double phi1 = A * cosh(k * (y + h)) / cosh(k * h) * sin(k * x);
  double phi2 = 3. * ak * A / (8. * alpa) * (sq(alpa) - 1.) * (sq(alpa) - 1.) *
                cosh(2.0 * k * (y + h)) * sin(2.0 * k * x) / cosh(2.0 * k * h);
  double phi3 = 1. / 64. * (sq(alpa) - 1.) * (sq(alpa) + 3.) *
                (9. * sq(alpa) - 13.) * cosh(3. * k * (y + h)) /
                cosh(3. * k * h) * a * a * k * k * A * sin(3. * k * x);
  Phi = phi1 + ak * phi2 + ak * ak * phi3;

  *p = Phi;
}

int main() {
  double x, y, vx, vy, eta, phi;

  x = 0.03;
  y = 0.04;

  f0(x, y, &vx, &vy, &eta);
  printf("%g %g %g %g %g\n", x, y, vx, vy, eta);

  f1(x, y, &vx, &vy, &eta);
  printf("%g %g %g %g %g\n", x, y, vx, vy, eta);

  f2(x, y, &phi);
  printf("%g %g %g\n", x, y, phi);
}
