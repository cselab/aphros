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
  double m[6];
  double v[3];

  argv = argv0;
  argv++;

  scl(&m[XX]);
  scl(&m[XY]);
  scl(&m[XZ]);
  scl(&m[YY]);
  scl(&m[YZ]);
  scl(&m[ZZ]);

  math_eig_values(m, v);

  printf("%g %g %g\n", v[X], v[Y], v[Z]);
  return 0;
}
