#include <vof.h>
#define MAX_SIZE (1024)
typedef struct Shape Shape;
struct Shape {
  double x, y, z;
  double rx, ry, rz;
  double u, v, w;
  int type;
};
enum {SHAPE_S, SHAPE_C};
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
  double d, p;
  x -= b.x;
  y -= b.y;
  z -= b.z;
  d = x*b.u + y*b.v + z*b.w;
  p = sq(b.u) + sq(b.v) + sq(b.w);
  if (!(p > 0)) {
    fprintf(stderr, "%s:%d: !(p > 0)\n", __FILE__, __LINE__);
    exit(2);
  }
  d /= sqrt(p);
  d = fabs(d);
  return d < b.rx;
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
