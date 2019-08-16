#include "fractions.h"
#include "curvature.h"

#include <unistd.h>

#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef POINTXYZ
#define POINTXYZ
#endif

#ifndef foreach_end
#define foreach_end
#endif

#ifndef foreach_neighbor_end
#define foreach_neighbor_end
#endif

#ifndef CELLIDX
#define CELLIDX
#endif

typedef struct {
  int Np;       // number of particles per string
  int Ns;       // number of strings (cross sections) per cell
  double Hp;    // lenght of string relative to cell size
  double eps;   // threshold for convergence
  int itermax;  // maximum number of iterations
  double eta;   // relaxation factor
  bool csv;     // dump csv particles and attraction points
  bool dumpit;  // dump the number of iterations and difference
  bool circ;    // replace line segments with circular arcs
} Partstr;

const int kMaxSection = 125 * 2; // maximum number of endpoints
                                 // from neighbor cells
                                 // (5x5x5 stencil and 2 endpoints in each)
const int kMaxFacet = 12;        // maximum number of vertices per facet

const int kMaxNp = 31;           // maximum value of kPartstr.Np


#if dimension == 2
#define kNs 1
#else 
#define kNs 3
#endif

static Partstr kPartstr = {7, kNs, 4., 1e-5, 20, 0.5, false, false, false};

#undef kNs

// Writes legacy vtk polydata
// fn: path
// xx: points
// nx: size of xx
// pp: polygons as lists of indices
// np: number of polygons
// ss[i] is size of polygon pp[i]
// cm: comment
// poly: true: polygons, false: lines
void WriteVtkPoly(const char* fn, coord* xx, int nx,
                  int* pp, int np, int* ss, const char* cm, bool poly) {
  FILE* o = fopen(fn, "w");
  fprintf(o, "# vtk DataFile Version 2.0\n");
  fprintf(o, "%s\n", cm);

  fprintf(o, "ASCII\n");
  fprintf(o, "DATASET POLYDATA\n");

  fprintf(o, "POINTS %d float\n", nx);
  for (int i = 0; i < nx; ++i) {
    coord x = xx[i];
    fprintf(o, "%g %g %g\n", x.x, x.y, x.z);
  }

  int na = 0; // total number of vortices
  for (int i = 0; i < np; ++i) {
    na += ss[i];
  }
  fprintf(o, "%s %d %d\n", poly ? "POLYGONS" : "LINES", np, np + na);
  int k = 0;
  for (int i = 0; i < np; ++i) {
    fprintf(o, "%d", ss[i]);
    for (int j = 0; j < ss[i]; ++j) {
      fprintf(o, " %d", pp[k++]);
    }
    fprintf(o, "\n");
  }

  fclose(o);
}

// Unit vector at angle ph.
static coord E(double ph) {
  coord p;
  p.x = cos(ph);
  p.y = sin(ph);
  p.z = 0;
  return p;
}


// Third component of cross product.
static double Cross3(coord a, coord b) {
  return a.x * b.y - a.y * b.x;
}

static coord Cross(coord a, coord b) {
  coord r;
  r.x = a.y*b.z - a.z*b.y;
  r.y = a.z*b.x - a.x*b.z;
  r.z = a.x*b.y - a.y*b.x;
  return r;
}

static coord Add(coord a, coord b) {
  a.x += b.x;
  a.y += b.y;
  a.z += b.z;
  return a;
}

static coord Sub(coord a, coord b) {
  a.x -= b.x;
  a.y -= b.y;
  a.z -= b.z;
  return a;
}

static coord Mul(coord a, double k) {
  a.x *= k;
  a.y *= k;
  a.z *= k;
  return a;
}

static coord Div(coord a, double k) {
  a.x /= k;
  a.y /= k;
  a.z /= k;
  return a;
}

static coord Coord(double x, double y, double z) {
  coord r = {x, y, z};
  return r;
}

static coord Zero() {
  coord r;
  r.x = 0;
  r.y = 0;
  r.z = 0;
  return r;
}

static double Dot(coord a, coord b) {
  double s = 0;
  s += a.x * b.x;
  s += a.y * b.y;
  s += a.z * b.z;
  return s;
}

static double Sqnorm(coord a) {
  return Dot(a, a);
}

static double Norm(coord a) {
  return sqrt(Sqnorm(a));
}

static double NormMax(coord a) {
  double r = fabs(a.x);
  if (fabs(a.y) > r) { r = fabs(a.y); }
  if (fabs(a.z) > r) { r = fabs(a.z); }
  return r;
}

static double Sqdist(coord a, coord b) {
  return Sqnorm(Sub(a, b));
}

static double Dist(coord a, coord b) {
  return sqrt(Sqdist(a, b));
}

static double Dotv(int np,
                   const coord* aa, const coord* bb) {
  double s = 0;
  for (int i = 0; i < np; ++i) {
    s += Dot(aa[i], bb[i]);
  }
  return s;
}

// Positions of particles.
// p: central particle
// ph: orientation angle (phi)
// th: bending angel (theta)
// np: number of particles
// hp: distance between particles
// Output:
// xx: array of length np
static void X(coord p, double ph, double th, int np, double hp, coord* xx) {
  int c = np / 2;
  xx[c] = p;
  for (int j = 0; j < c; ++j) {
    double jp = j + 0.5;
    xx[c + j + 1] = Add(xx[c + j], Mul(E(ph + jp * th), hp));
    xx[c - j - 1] = Add(xx[c - j], Mul(E(ph - jp * th), -hp));
  }
}

// Point nearest to line segment [a, b]
static coord Nearest(coord a, coord b, coord x) {
  b = Sub(b, a);
  x = Sub(x, a);
  double q = Dot(b, x) / Dot(b, b);
  q = clamp(q, 0., 1.);
  a.x += b.x * q;
  a.y += b.y * q;
  a.z = 0;
  return a;
}

// Derivative dX/dph
static void DxDph(coord p, double ph, double th, 
                  int np, double hp, coord* xx) {
  (void) p;
  int c = np / 2;
  xx[c] = Zero();
  for (int j = 0; j < c; ++j) {
    double jp = j + 0.5;
    xx[c + j + 1] = Add(xx[c + j],
                        Mul(E(ph + jp * th + PI * 0.5), hp));
    xx[c - j - 1] = Add(xx[c - j],
                        Mul(E(ph - jp * th + PI * 0.5), -hp));
  }
}

// Derivative dX/dth
static void DxDth(coord p, double ph, double th, 
                  int np, double hp, coord* xx) {
  (void) p;
  int c = np / 2;
  xx[c] = Zero();
  for (int j = 0; j < c; ++j) {
    double jp = j + 0.5;
    xx[c + j + 1] = Add(xx[c + j],
                        Mul(E(ph + jp * th + PI * 0.5), hp * jp));
    xx[c - j - 1] = Add(xx[c - j],
                        Mul(E(ph - jp * th + PI * 0.5), hp * jp));
  }
}

// Displacement from line to circular arc.
// k: curvature
// w: distance from endpoint to center
// d: distance from point to center
static double CircDelta(double k, double w, double d) {
  const double th = 0.99;
  k = clamp(k, -th / w, th / w); // limit the curvature for sqrt and division
                                 // (w must be smaller than radius of curvaure)
  d = min(d, w);
  double t1 = sqrt(1. - sq(k * w));
  double t2 = sqrt(1. - sq(k * d));
  return k * (sq(w) - sq(d)) / (t1 + t2);
}

// Point on circular arc.
// k: curvature of circular arc
// a, b: endpoints
// y: point on line
static coord Circ(double k, coord a, coord b, coord y) {
  coord q = Sub(b, a);
  coord n;
  n.x = q.y / Norm(q);
  n.y = -q.x / Norm(q);
  n.z = 0;
  coord c = Mul(Add(a, b), 0.5);
  double d = Dist(c, y);
  double w = Dist(c, a);
  return Add(y, Mul(n, CircDelta(k, w, d)));
}

// Forces on particles.
// np: number of particles
// xx: positions 
// nl: number of points in ll
// ll: flat array of endpoints of line segments
// eta: relaxation factor
// k: curvature
// Output:
// ff: forces
static void F(int np, const coord* xx, int nl, const coord* ll, 
              double eta, double k, bool circ,coord* ff) {
  if (nl == 0) {
    for (int i = 0; i < np; ++i) {
      ff[i] = Zero();
    }
    return;
  }
  for (int i = 0; i < np; ++i) {
    int lm = 0;
    coord pm = Nearest(ll[0], ll[1], xx[i]);

    for (int l = 0; l < nl; l += 2) {
      coord p = Nearest(ll[l], ll[l + 1], xx[i]);
      if (Sqdist(xx[i], p) < Sqdist(xx[i], pm)) {
        lm = l;
        pm = p;
      }
    }
    if (circ) {
      pm = Circ(k, ll[lm], ll[lm + 1], pm);
    }
    ff[i] = Mul(Sub(pm, xx[i]), eta);
  }
}

// Iteration of evolution of particles.
// p_,ph_,th_: current configuration
// np: number of particles
// hp: distance between particles
// ff: forces (output of F())
// xx: buffer of size np
// Output:
// p_,ph_,th_: new configuration
// ff, xx: modified
// Returns maximum absolute difference.
static double Iter(coord* p_, double* ph_, double* th_,
                 int np, double hp, coord* ff, coord* xx) {
  coord p = *p_;
  double ph = *ph_;
  double th = *th_;

  coord tt[kMaxNp];
  int c = np / 2;

  X(p, ph, th, np, hp, xx);

  // correct p
  p = Add(p, ff[c]);
  X(p, ph, th, np, hp, tt);
  for (int i = 0; i < np; ++i) {
    coord dx = Sub(tt[i], xx[i]);
    xx[i] = Add(xx[i], dx);
    ff[i] = Sub(ff[i], dx);
  }

  // correct phi
  DxDph(p, ph, th, np, hp, tt);
  ph += Dotv(np, tt, ff) / Dotv(np, tt, tt);
  X(p, ph, th, np, hp, tt);
  for (int i = 0; i < np; ++i) {
    coord dx = Sub(tt[i], xx[i]);
    xx[i] = Add(xx[i], dx);
    ff[i] = Sub(ff[i], dx);
  }

  // correct theta
  DxDth(p, ph, th, np, hp, tt);
  th += Dotv(np, tt, ff) / Dotv(np, tt, tt);
  X(p, ph, th, np, hp, tt);

  double r = 0;
  for (int i = 0; i < np; ++i) {
    r = max(r, NormMax(Sub(xx[i], tt[i])));
  }

  *p_ = p;
  *ph_ = ph;
  *th_ = th;
  return r;
}

// incid: if true increment id
// k: attribute
static void AppendCsv(int it, int np, coord * xx, const char* name, 
                      double c, double k) {
  FILE* o;
  char on[255];
  sprintf(on, "%s_%04d.csv", name, it);
  if (access(on, W_OK) == -1) { // new file
    o = fopen(on, "w");
    fprintf(o, "x,y,z,c,k\n");
  } else { // append
    o = fopen(on, "a");
  }

  for (int i = 0; i < np; ++i) {
    fprintf(o, "%g,%g,%g,%g,%g", xx[i].x, xx[i].y, xx[i].z, c, k);
    fprintf(o, "\n");
  }
  fclose(o);
}

static double Curv(double hp, double th) {
  return sqrt(2.) / hp * sin(th) / sqrt(1. + cos(th));
}

static void Swap(coord* a, coord* b) {
  coord t;
  t = *a;
  *a = *b;
  *b = t;
}


static coord Abs(coord p) {
  p.x = fabs(p.x);
  p.y = fabs(p.y);
  p.z = fabs(p.z);
  return p;
}

static double Get(coord p, int i) {
  return i == 0 ? p.x : i == 1 ? p.y : p.z;
}

// Index of minimal component.
static int Argmin(coord p) {
  int im = 0;
  for (int i = 1; i < 3; ++i) {
    if (Get(p, i) < Get(p, im)) {
      im = i;
    }
  }
  return im;
}

// Index of maximum component.
static int Argmax(coord p) {
  int im = 0;
  for (int i = 1; i < 3; ++i) {
    if (Get(p, i) > Get(p, im)) {
      im = i;
    }
  }
  return im;
}

static void Set(coord* p, int i, double a) {
  if (i == 0) { p->x = a; }
  else if (i == 1) { p->y = a; }
  else { p->z = a; }
}

// Returns 3D base aligned with mesh directions.
// n: unit vector
// Output:
// *t, *u: vectors to form orthonormal base <t,u,n>
static void GetBase(coord n, coord* t, coord* u) {
  int i = Argmin(Abs(n));
  coord e = Zero();
  Set(&e, i, 1.);
  n = Div(n, Norm(n));
  *t = Cross(n, e);
  *t = Div(*t, Norm(*t));
  *u = Cross(*t, n);
}

// Transformation of coordinates.
//    |n
//    |
//    |
//    o-----t
//   /
//  /u
typedef struct {
  coord o; // origin
  coord t, n, u; // orthonormal positive-oriented base
} Trans;

// p: point
// o: origin
// t,n,u: orthonormal base
// Output:
// *l: coordinates of p in <o,t,n,u>
static coord GlbToLoc(coord p, Trans w) {
  p = Sub(p, w.o);
  coord l = {Dot(w.t, p), Dot(w.n, p), Dot(w.u, p)};
  return l;
}

// l: coordinates in <t,n,u,o>
// t,n,u: orthonormal base
// o: origin
// Output:
// *p: point
static coord LocToGlb(coord l, Trans w) {
  coord p = w.o;
  p = Add(p, Mul(w.t, l.x));
  p = Add(p, Mul(w.n, l.y));
  p = Add(p, Mul(w.u, l.z));
  return p;
}

// Returns the number of interfacial cells.
static int GetNcInter(scalar c) {
  int nc = 0; 
  foreach() {
    if (interfacial(point, c)) {
      ++nc;
    }
  }
  foreach_end
  return nc;
}

// Set z-copmonent to zero if 2D.
static void A2(coord* p) {
  (void) p;
#if dimension == 2
  p->z = 0;
#endif
}

static int Facets(coord m, double alpha, coord* pp) {
#if dimension == 2
  int nf = facets(m, alpha, pp);
  for (int i = 0; i < nf; ++i) {
    A2(&pp[i]);
  }
  return nf;
#else
  return facets(m, alpha, pp, 1.);
#endif
}

static coord Mycs(Point point, scalar c) {
  coord m = mycs(point, c);
  A2(&m);
  return m;
}

// Cross-section of 2D interface from neighbor cells in plane coordinates.
// point: center cell
// c: volume fraction
// a: transformation
// ll: buffer for at least kMaxSection more points
// *nl: current size of ll
// Output:
// ll: appended by local coordinates of endpoints,
//     p = o + t*l.x  + t*l.y ,  [pa,pb] is one line segment
// *nl: new size of ll
static void Section2(Point point, scalar c, vector nn, Trans w,
                    coord* ll, int* nl) {
  foreach_neighbor(2)
  {
    if (c[CELLIDX] > 0. && c[CELLIDX] < 1.) {
      coord m = {nn.x[CELLIDX], nn.y[CELLIDX], nn.z[CELLIDX]};
      double alpha = plane_alpha(c[CELLIDX], m);
      coord pp[kMaxFacet];
      int nf = Facets(m, alpha, pp);
      coord rn = {x,y,z};
      if (nf == 2) {
        coord p = Add(rn, Mul(pp[0], Delta));
        coord pb = Add(rn, Mul(pp[1], Delta));
        coord l = GlbToLoc(p, w);
        coord lb = GlbToLoc(pb, w);

        coord dl = Sub(l, lb);
        double mt = Dot(w.t, m);
        double mn = Dot(w.n, m);
        double dt = dl.x;
        double dn = dl.y;
        if (Cross3(Coord(mt, mn, 0), Coord(dt, dn, 0)) < 0) {
          Swap(&l, &lb);
        }
        if (Sqnorm(dl) > 0) {
          ll[*nl] = l;
          ll[*nl + 1] = lb;
          *nl += 2;
        }
      }
    }
  }
  foreach_neighbor_end
}

// Cross-section of 3D interface from neighbor cells in plane coordinates.
// point: center cell
// c: volume fraction
// w: transformation
// ll: buffer for at least kMaxSection more points
// *nl: current size of ll
// Output:
// ll: appended by local coordinates of endpoints,
//     p = o + t*l.x  + t*l.y ,  [pa,pb] is one line segment
// *nl: new size of ll
static void Section3(Point point, scalar c, vector nn, Trans w,
                    /**/ coord* ll, int* nl) {
  foreach_neighbor(2)
  {
    if (c[CELLIDX] > 0. && c[CELLIDX] < 1.) {
      coord m = {nn.x[CELLIDX], nn.y[CELLIDX], nn.z[CELLIDX]};
      double alpha = plane_alpha(c[CELLIDX], m);
      coord pp[kMaxFacet];
      int nf = Facets(m, alpha, pp);
      coord rn = {x,y,z};
      assert(nf <= kMaxFacet);
      int q = 0;  // number of intersections found

      for (int i = 0; i < nf && q < 2; ++i) {
        int ib = (i + 1) % nf;
        coord p = Add(rn, Mul(pp[i], Delta));
        coord pb = Add(rn, Mul(pp[ib], Delta));
        coord l = GlbToLoc(p, w);
        coord lb = GlbToLoc(pb, w);

        // intersection with l.z=0:
        //   0 = l.z * s + lb.z * (1 - s)
        if (l.z * lb.z < 0) { // line crosses the plane
          double s = lb.z / (lb.z - l.z);
          ll[*nl + q].x = l.x * s + lb.x * (1 - s);
          ll[*nl + q].y = l.y * s + lb.y * (1 - s);
          ll[*nl + q].z = 0.;
          if (++q == 2) {
            coord dl = Sub(ll[*nl + 1], ll[*nl]);
            double mt = Dot(w.t, m);
            double mn = Dot(w.n, m);
            double dt = dl.x;
            double dn = dl.y;
            if (Cross3(Coord(mt, mn, 0), Coord(dt, dn, 0)) > 0) {
              Swap(&ll[*nl], &ll[*nl + 1]);
            }
            if (Sqnorm(dl) > 0) {
              *nl += 2;
            }
          }
        }
      }
    }
  }
  foreach_neighbor_end
}


static void Section(Point point, scalar c, vector nn,
                    Trans w, coord* ll, int* nl) {
#if dimension == 2
  Section2(point, c, nn, w, ll, nl);
#else
  Section3(point, c, nn, w, ll, nl);
#endif
}

// Curvature of a set line segments.
// ll: flat array endpoints of line segments
// nl: size of ll
// *w: transformation (for csv output)
// res_: difference at last iteration
// it_: number of iterations
// Output:
// appends positions and attraction points to csv if a is not null
static double GetLinesCurv(coord* ll, int nl, double delta, const Trans* w,
                           Partstr conf, double* res_, int* it_, 
                           double hash) {
  if (nl >= 4) { // require at least two segments
    const int Np = conf.Np;
    const double eta = conf.eta;
    const int itermax = conf.itermax;
    const double hp = conf.Hp * delta / (Np - 1);

    coord xx[kMaxNp];  // positions
    coord ff[kMaxNp];  // forces
    coord p = {0., 0, 0.};
    double ph = 0.;
    double th = 0.;
    double k = 0;
    double res = 0;
    int it = 0;
    for (it = 0; it < itermax; ++it) {
      X(p, ph, th, Np, hp, xx);
      F(Np, xx, nl, ll, eta, k, conf.circ, ff);

      k = Curv(hp, th);
      res = Iter(&p, &ph, &th, Np, hp, ff, xx);
      if (res / (eta * delta) < conf.eps) {
        break;
      }
    }

    if (w) {
      coord tt[kMaxNp];
      for (int i = 0; i < Np; ++i) {
        tt[i] = Add(xx[i], Mul(ff[i], 1. / eta));
        tt[i] = LocToGlb(tt[i], *w);
      }
      for (int i = 0; i < Np; ++i) {
        xx[i] = LocToGlb(xx[i], *w);
      }

      AppendCsv(t * 100, Np, xx, "part", hash, -k);
      AppendCsv(t * 100, Np, tt, "attr", hash, -k);
    }

    *res_ = res;
    *it_ = it;
    if (conf.dumpit) {
      fprintf(stderr, "%d %g \n", it, res);
    }
    return k;
  }
  return 0;
}

// Curvature in cross section by plane through a.n,a.t and point a.o
// point: target cell
static double GetCrossCurv(Point point, scalar c, vector nn,
                           Trans w, Partstr conf) {
  coord ll[kMaxSection];
  int nl = 0;
  Section(point, c, nn, w, ll, &nl);
  double res;
  int it;
#ifndef NOBA
#if dimension == 2
  double hash = 1000 * point.j + point.i;
#else
  double hash = 1000 * (1000 * point.k + point.j) + point.i;
#endif
#else
  static double hash = 0;
  hash += 1.;
#endif
  return GetLinesCurv(
      ll, nl, Delta, conf.csv ? &w : NULL, conf, &res, &it, hash);
}

// Transformation b rotated at angle.
// s: index of cross sections
// Ns: number of cross sections
static Trans GetSectionTrans(int s, int Ns, Trans b) {
  const double g = PI *  s / Ns;
  Trans w = b;
  w.t = Add(Mul(b.t, cos(g)), Mul(b.u, sin(g)));
  w.u = Cross(w.t, w.n);
  return w;
}

// Mean curvature over multiple cross sections by planes rotates around b.n
static double GetMeanCurv(Point point, scalar c, vector nn,
                          Trans b, Partstr conf) {
  double ksum = 0;
  const int Ns = conf.Ns;
  for (int s = 0; s < Ns; ++ s) {
    Trans w = GetSectionTrans(s, Ns, b);
    ksum += GetCrossCurv(point, c, nn, w, conf);
  }
  return ksum / Ns;
}

// Dumps interface fragments to vtk.
void DumpFacets(scalar c, const char* fn) {
  const int nc = GetNcInter(c); // number of interfacial cells
  const int mm = nc * kMaxFacet;

  coord* xx = (coord*)malloc(mm * sizeof(coord));
  int* pp = (int*)malloc(mm * sizeof(int));
  int* ss = (int*)malloc(mm * sizeof(int));
  int np = 0;
  int nx = 0;

  foreach() {
    if (interfacial (point, c)) {
      coord m = Mycs(point, c);
      double alpha = plane_alpha (c[CELLIDX], m);
      coord p[kMaxFacet];
      int nf = Facets(m, alpha, p);

      coord r = {x, y, z};

      ss[np] = nf;
      for (int i = 0; i < nf; ++i) {
        pp[nx] = nx;
        xx[nx] = Add(r, Mul(p[i], Delta));
        ++nx;
      }
      ++np;
    }
  }
  foreach_end
  WriteVtkPoly(fn, xx, nx, pp, np, ss, "comment", true);

  free(ss);
  free(pp);
  free(xx);
}

// Returns transformation to local coordinates at the interface.
static Trans GetPointTrans(Point point, scalar c, vector nn) {
  Trans b;
  coord n = {nn.x[CELLIDX], nn.y[CELLIDX], nn.z[CELLIDX]};
  double alpha = plane_alpha(c[CELLIDX], n);

  coord o;
  plane_area_center(n, alpha, &o);
  A2(&o);
  POINTXYZ
  coord rc = {x, y, z};
  o = Add(rc, Mul(o, Delta));

  b.o = o;
  b.n = Div(n, Norm(n));
  GetBase(b.n, &b.t, &b.u);
  return b;
}

static void CalcNormal(scalar c, vector nn) {
  foreach() {
    coord n = Mycs(point, c);
    nn.x[CELLIDX] = n.x;
    nn.y[CELLIDX] = n.y;
    nn.z[CELLIDX] = n.z;
  }
  foreach_end
}

// Dumps cross sections of the inteface fragments to vtk.
void DumpLines(scalar c, vector nn, Partstr conf, const char* fn) {
  const int Ns = conf.Ns;

  const int nc = GetNcInter(c); // number of interfacial cells
  const int mm = nc * kMaxSection;

  coord* xx = (coord*)malloc(mm * sizeof(coord));
  int* pp = (int*)malloc(mm * sizeof(int));
  int* ss = (int*)malloc(mm * sizeof(int));
  int nx = 0;
  int np = 0;

  foreach() {
    if (interfacial(point, c)) {
      Trans b = GetPointTrans(point, c, nn);

      for (int s = 0; s < Ns; ++ s) {
        Trans w = GetSectionTrans(s, Ns, b);

        coord ll[kMaxSection];
        int nl = 0;

        Section(point, c, nn, w, ll, &nl);
        GetCrossCurv(point, c, nn, w, conf);

        for (int i = 0; i < nl; ++i) {
          ll[i] = LocToGlb(ll[i], w);
        }

        for (int i = 0; i < nl; i += 2) {
          ss[np] = 2;
          pp[nx] = nx;
          xx[nx] = ll[i];
          ++nx;
          pp[nx] = nx;
          xx[nx] = ll[i + 1];
          ++nx;
          ++np;
        }
      }
    }
  }
  foreach_end
  WriteVtkPoly(fn, xx, nx, pp, np, ss, "lines", false);

  free(ss);
  free(pp);
  free(xx);
}

static double partstr(Point point, scalar c, vector nn) {
  Trans b = GetPointTrans(point, c, nn);
  double k = -GetMeanCurv(point, c, nn, b, kPartstr);
#if dimension == 3
  k *= 2;
#endif
  return k;
}

#ifndef NOBA
trace
cstats curvature_partstr(struct Curvature p)
{
  scalar c = p.c, kappa = p.kappa;
  double sigma = p.sigma ? p.sigma : 1.;
  int sh = 0, sf = 0, sa = 0, sc = 0;
  vector ch = c.height, h = automatic (ch);
  if (!ch.x.i)
    heights (c, h);

#if TREE
  kappa.refine = kappa.prolongation = curvature_prolongation;
  kappa.restriction = curvature_restriction;
#endif

  scalar k[];
  scalar_clone (k, kappa);

  vector nn[];
  CalcNormal(c, nn);

  boundary({nn});

  foreach(reduction(+:sh) reduction(+:sc)) {
    if (!interfacial (point, c)) {
      k[CELLIDX] = nodata;
#ifndef PS_nohf
    } else if ((k[CELLIDX] = height_curvature (point, c, h)) != nodata) {
      sh++;
#endif
    } else  {
      k[CELLIDX] = partstr(point, c, nn);
      sc++;
    }
  }
  foreach_end
  boundary ({k});

  foreach () {
    double kf = k[CELLIDX];
    if (kf == nodata)
      kappa[CELLIDX] = nodata;
    else if (p.add)
      kappa[CELLIDX] += sigma*kf;
    else
      kappa[CELLIDX] = sigma*kf;
  }
  foreach_end
  boundary ({kappa});

  return (cstats){sh, sf, sa, sc};
}

trace
cstats curvature_orig(struct Curvature p) {
  return curvature(p);
}
#define curvature curvature_partstr

#endif

