#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inside.h>

static const char* me = "plane";
static double eps = 1e-6;
#define USED(x) \
  if (x)        \
    ;           \
  else {        \
  }

static void usg(void) {
  fprintf(stderr, "%s -a alpha -n x y z mesh | ./2off > off\n", me);
  exit(1);
}
static int plain_tri(
    const double[3], double, const double[3], const double[3], const double[3],
    double*);
static double dist2(const double[3], const double[3]);
int main(int argc, const char** argv) {
  USED(argc);
  enum { X, Y, Z };
  const double* a;
  const double* b;
  const double* c;
  double alpha;
  double n[3];
  double p[9];
  double* pp;
  double* ver;
  int Aflag;
  int cap;
  int cnt;
  int Found;
  int i;
  int j;
  int k;
  int l;
  int Nflag;
  int np;
  int nt;
  int nv;
  int face[3];
  int t;
  int* tri;

  Aflag = Nflag = 0;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
      case 'h':
        usg();
        break;
      case 'a':
        argv++;
        if (*argv == NULL) {
          fprintf(stderr, "%s: -a needs an argument\n", me);
          return 2;
        }
        alpha = atof(*argv);
        Aflag = 1;
        break;
      case 'n':
        argv++;
        if (argv[0] == NULL || argv[1] == NULL || argv[2] == NULL) {
          fprintf(stderr, "%s: -n needs three arguments\n", me);
          return 2;
        }
        n[X] = atof(*argv++);
        n[Y] = atof(*argv++);
        n[Z] = atof(*argv);
        if (n[X] == 0 && n[Y] == 0 && n[Z] == 0) {
          fprintf(stderr, "%s: normal cannot be zero\n", me);
          return 2;
        }
        Nflag = 1;
        break;
      default:
        fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
        exit(2);
    }
  if (argv[0] == NULL) {
    fprintf(stderr, "%s: mesh file is missing\n", me);
    exit(2);
  }
  if (Nflag == 0) {
    fprintf(stderr, "%s: -n is not set\n", me);
    exit(2);
  }
  if (Aflag == 0) {
    fprintf(stderr, "%s: -a is not set\n", me);
    exit(2);
  }

  if (inside_mesh_read(argv[0], &nt, &tri, &nv, &ver) != 0) {
    fprintf(stderr, "%s: fail to read mesh '%s'\n", me, argv[0]);
    exit(2);
  }
  np = cap = 0;
  pp = NULL;
  for (t = 0; t < nt; t++) {
    i = tri[3 * t];
    j = tri[3 * t + 1];
    k = tri[3 * t + 2];
    a = &ver[3 * i];
    b = &ver[3 * j];
    c = &ver[3 * k];
    cnt = plain_tri(n, alpha, a, b, c, p);
    if (cnt) {
      for (k = 0; k < cnt; k++) {
        Found = 0;
        for (l = 0; l < np; l++)
          if (dist2(&p[3 * k], &pp[3 * l]) < eps * eps) {
            Found = 1;
            break;
          }
        if (Found) {
        } else {
          if (np >= cap) {
            cap = 2 * cap + 1;
            if ((pp = realloc(pp, 3 * cap * sizeof *pp)) == NULL) {
              fprintf(stderr, "%s: realloc failed\n", me);
              return 2;
            }
          }
          l = np;
          pp[3 * np] = p[3 * k];
          pp[3 * np + 1] = p[3 * k + 1];
          pp[3 * np + 2] = p[3 * k + 2];
          np++;
        }
        face[k] = l;
      }
      printf("f");
      for (k = 0; k < cnt; k++) {
        printf(" %d", face[k]);
      }
      printf("\n");
    }
  }

  for (i = 0; i < np; i++)
    printf("p %-.16e %-.16e %-.16e\n", pp[3 * i], pp[3 * i + 1], pp[3 * i + 2]);

  free(pp);
  inside_mesh_fin(tri, ver);
}

static double plane_point(const double n[3], double alpha, const double p[3]) {
  enum { X, Y, Z };
  return n[X] * p[X] + n[Y] * p[Y] + n[Z] * p[Z] + alpha;
}

static int plain_edg(
    const double n[3], double alpha, const double a[3], const double b[3],
    /**/ double p[3]) {
  enum { X, Y, Z };
  double x;
  double y;
  double t;
  x = plane_point(n, alpha, a);
  if (fabs(x) < eps) {
    p[X] = a[X];
    p[Y] = a[Y];
    p[Z] = a[Z];
    return 1;
  }
  y = plane_point(n, alpha, b);
  if (fabs(y) < eps) {
    p[X] = b[X];
    p[Y] = b[Y];
    p[Z] = b[Z];
    return 1;
  }
  if (x * y > 0) return 0;
  t = x / (x - y);
  p[X] = a[X] + t * (b[X] - a[X]);
  p[Y] = a[Y] + t * (b[Y] - a[Y]);
  p[Z] = a[Z] + t * (b[Z] - a[Z]);
  return 1;
}

static int plain_tri(
    const double n[3], double alpha, const double a[3], const double b[3],
    const double c[3], double* ver) {
  int cnt;
  cnt = 0;
  if (plain_edg(n, alpha, a, b, ver)) {
    cnt++;
    ver += 3;
  }
  if (plain_edg(n, alpha, b, c, ver)) {
    cnt++;
    ver += 3;
  }
  if (plain_edg(n, alpha, c, a, ver)) {
    cnt++;
    ver += 3;
  }
  return cnt;
}

static double sq(double x) {
  return x * x;
}
static double dist2(const double a[3], const double b[3]) {
  enum { X, Y, Z };
  return sq(a[X] - b[X]) + sq(a[Y] - b[Y]) + sq(a[Z] - b[Z]);
}
