#include <stdio.h>
#include <stdlib.h>
#include <inside.h>
#include <math.h>

const char* me = "circumradius";

static void
usg(void) {
  fprintf(stderr, "%s -i float float float float float float float float float\n", me);
  exit(1);
}

static double
sq(double x)
{
  return x*x;
}
static double
edg(const double *a, const double *b)
{
  enum {X, Y, Z};
  return sqrt(sq(a[X] - b[X]) + sq(a[Y] - b[Y]) + sq(a[Y] - b[Y]));
}
static int
circumradius(const double * u, const double * v, const double * w, double *r)
{
  double a;
  double b;
  double c;
  double s;
  double num;
  double den;
  a = edg(v, u);
  b = edg(w, u);
  c = edg(v, w);
  s = (a + b + c)/2;
  num = a*b*c;
  den = s*(s - a)*(s - b)*(s - c);
  if (den <= 0)
    return 1;
  *r = num/(4*sqrt(den));
  return 0;
}

int
main(int argc, const char** argv) {
  enum {X, Y, Z};
  double a[3];
  double b[3];
  double c[3];
  double r;

  argv++; argc--;
  while (*argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
    case 'h':
      usg();
      break;
    case 'i':
      argv++; argc--;
      if (argc < 9) {
	fprintf(stderr, "%s: needs nine arguments given '%d'\n", me, argc);
	exit(2);
      }
      a[X] = atof(*argv++); argc--;
      a[Y] = atof(*argv++); argc--;
      a[Z] = atof(*argv++); argc--;
      b[X] = atof(*argv++); argc--;
      b[Y] = atof(*argv++); argc--;
      b[Z] = atof(*argv++); argc--;      
      c[X] = atof(*argv++); argc--;
      c[Y] = atof(*argv++); argc--;
      c[Z] = atof(*argv++); argc--;      
      break;
    default:
      fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
      exit(2);
    }
  if (argv[0] != NULL) {
      fprintf(stderr, "%s: unknown argument '%s'\n", me, argv[0]);
      exit(2);    
  }
  if (circumradius(a, b, c, &r) != 0) {
    fprintf(stderr, "%s: circumradius failed\n", me);
    exit(2);
  }
  printf("%.16g\n", r);
}
