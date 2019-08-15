#include "grid/multigrid3D.h"
#include "fractions.h"
#include "curvature_select.h"
#include "io/io.h"
#include <vof.h>

#include <assert.h>
#define myassert(EX) (void)((EX) || (__assert (#EX, __FILE__, __LINE__),0))

#include <mpi.h>

#define ONROOT int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); if (rank == 0)
#define MAX_SIZE (1024)

typedef struct Shape Shape;
struct Shape {
  double x, y, z;
  double rx, ry, rz;
  double u, v, w;
  int type;
};
enum {SHAPE_S, SHAPE_C};

int nxexp;
int argnx;
void CreateFieldS(Shape b, scalar c) {
  foreach() {
    double h = Delta;
    c[] = vof_3d((b.x - x)/h , (b.y - y)/h, (b.z - z)/h, b.rx/h);
  }
}
int GoodS(Shape b, double x, double y, double z, double h) {
  return 1;
}
void CreateFieldC(Shape b, scalar c) {
  foreach() {
    double h = Delta;
    c[] = vof_cylinder((b.x - x)/h , (b.y - y)/h, (b.z - z)/h,
                       b.rx/h, b.u, b.v, b.w);
  }
}
int GoodC(Shape b, double x, double y, double z, double h) {
  return sq(x -  b.x) + sq(y -  b.y) + sq(z -  b.z) < sq(b.rx + h);
}
void (*CreateFieldArray[])(Shape, scalar) = {CreateFieldS, CreateFieldC};
int (*GoodArray[])(Shape, double, double, double, double) = {GoodS, GoodC};
void CreateField(Shape b, scalar c) {
  return CreateFieldArray[b.type](b, c);
}
int Good(Shape b, double x, double y, double z, double h) {
  return GoodArray[b.type](b, x, y, z, h);
}

Shape Read(FILE *f) {
  char t;
  Shape b;
  int n;
  char ss[MAX_SIZE], *s = ss;
  fgets(s, MAX_SIZE, f);
  sscanf(s, "%c", &t);
  if (('0' <= t && t <= '9') || t == '.' || t == '-' || t == '+')
    t = 's';
  else
    s++;
  switch (t) {
  case 's':
    b.type = SHAPE_S;
    n = sscanf(s, "%lf %lf %lf %lf", &b.x, &b.y, &b.z, &b.rx);
    if (n != 4) {
      fprintf(stderr, "%s:%d: wrong line '%s'\n", __FILE__, __LINE__, ss);
      exit(2);
    }      
    break;
  case 'c':
    b.type = SHAPE_C;
    n = sscanf(s, "%lf %lf %lf %lf %lf %lf %lf",
               &b.x, &b.y, &b.z, &b.rx, &b.u, &b.v, &b.w);
    if (n != 7) {
      fprintf(stderr, "%s:%d: wrong line '%s'\n", __FILE__, __LINE__, ss);
      exit(2);
    }
    break;
  default:
    fprintf(stderr, "%s:%d: unknown type '%c' in line '%s'\n",
            __FILE__, __LINE__, t, s);
    exit(2);
  }
  return b;
}

int main() {
  Shape b;
  {
    FILE* q = fopen("nxexp", "r");
    fscanf(q, "%d", &nxexp);
    fclose(q);
  }
  argnx = (1 << nxexp);
  init_grid(argnx);
  FILE* fb = fopen("b.dat", "r");
  b = Read(fb);
  fclose(fb);
  origin (0.,0.,0.);
  scalar vf[];
  CreateField(b, vf);
  scalar k[];
#ifdef CURV_PARTSTR
#ifdef PS_Np
  kPartstr.Np = PS_Np;
#endif
#ifdef PS_Ns
#if dimension == 3
  kPartstr.Ns = PS_Ns;
#endif
#endif
#ifdef PS_Hp
  kPartstr.Hp = PS_Hp;
#endif
#ifdef PS_eps
  kPartstr.eps = PS_eps;
#endif
#ifdef PS_itermax
  kPartstr.itermax = PS_itermax;
#endif
#ifdef PS_eta
  kPartstr.eta = PS_eta;
#endif
#ifdef PS_circ
  kPartstr.circ = PS_circ;
#endif
#endif
  curvature(vf, k);
#if dimension == 2
  double kc = 1.;
#elif dimension == 3
  double kc = 0.5;
#endif
  
  {
    FILE* q = fopen("ok", "w");
    foreach() {
      if (vf[] > 0. && vf[] < 1.) {
        if (Good(b, x, y, z, Delta))
          fprintf(q, "%.16g\n", k[] * kc);
      }
    }
    fclose(q);
  }

  {
    FILE* q = fopen("ovf", "w");
    double s = 0.;
    foreach() {
      double h = Delta;
      s += vf[] * h * h * h;
    }
    fprintf(q, "%.16g\n", s);
    fclose(q);
  }


  {
    FILE* q = fopen("u.vtk", "w");
    io({vf, k}, q);
    fclose(q);
  }
}
