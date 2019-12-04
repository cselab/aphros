#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <poc/math.h>

static const char *me = "vectors";

static const char **argv;
static int
scl( /**/ double *p)
{
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

int
main(int argc, const char **argv0)
{
  enum { X, Y, Z };
  enum { XX, XY, XZ, YY, YZ, ZZ };
  int i;

  double m[6] = { 1, 1 / 2.0, 1 / 3.0,
    1 / 3.0, 1 / 4.0,
    1 / 5.0
  };
  double v[3 * 3];

  argv = argv0;
  argv++;

  scl(&m[XX]);
  scl(&m[XY]);
  scl(&m[XZ]);
  scl(&m[YY]);
  scl(&m[YZ]);
  scl(&m[ZZ]);

  math_eig_vectors(m, v);

  i = 0;
  printf("%g %g %g\n", v[i + X], v[i + Y], v[i + Z]);
  i += 3;
  printf("%g %g %g\n", v[i + X], v[i + Y], v[i + Z]);
  i += 3;
  printf("%g %g %g\n", v[i + X], v[i + Y], v[i + Z]);
  i += 3;
  return 0;
}
