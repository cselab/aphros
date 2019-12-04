#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <poc/math.h>

static const char *me = "orto";
enum { X, Y, Z };
enum { XX, XY, XZ, YY, YZ, ZZ };
static double dot(double a[3], double b[3]);

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
  int i;
  double v[3 * 3], m[6];
  double a[3], b[3], c[3];

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
  a[X] = v[X];
  a[Y] = v[Y];
  a[Z] = v[Z];
  b[X] = v[3 + X];
  b[Y] = v[3 + Y];
  b[Z] = v[3 + Z];
  c[X] = v[6 + X];
  c[Y] = v[6 + Y];
  c[Z] = v[6 + Z];

  printf("%.3f %.3f %.3f\n", dot(a, a), dot(a, b), dot(a, c));
  printf("%.3f %.3f %.3f\n", dot(b, a), dot(b, b), dot(b, c));
  printf("%.3f %.3f %.3f\n", dot(c, a), dot(c, b), dot(c, c));

  return 0;
}

static double
dot(double a[3], double b[3])
{
  return a[X] * b[X] + a[Y] * b[Y] + a[Z] * b[Z];
}
