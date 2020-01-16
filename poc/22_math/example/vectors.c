// Created by Sergey Litvinov on 04.12.2019
// Copyright 2019 ETH Zurich

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <poc/math.h>

static const char* me = "vectors";

static const char** argv;
static int scl(/**/ double* p) {
  if (*argv == NULL) {
    fprintf(stderr, "%s: not enough args\n", me);
    exit(2);
  }
  if (sscanf(*argv, "%lf", p) != 1) {
    fprintf(stderr, "%s: not a number '%s'\n", me, *argv);
    exit(2);
  }
  argv++;
  return 0;
}

int main(int argc, const char** argv0) {
  enum { X, Y, Z };
  enum { XX, XY, XZ, YY, YZ, ZZ };
  int i;
  double a[3], b[3], c[3], m[6];

  argv = argv0;
  argv++;

  scl(&m[XX]);
  scl(&m[XY]);
  scl(&m[XZ]);
  scl(&m[YY]);
  scl(&m[YZ]);
  scl(&m[ZZ]);

  math_eig_vectors(m, a, b, c);

  printf("%g %g %g\n", a[X], a[Y], a[Z]);
  printf("%g %g %g\n", b[X], b[Y], b[Z]);
  printf("%g %g %g\n", c[X], c[Y], c[Z]);
  return 0;
}
