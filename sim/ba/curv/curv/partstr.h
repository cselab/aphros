#include "fractions.h"
#include "curvature.h"


#include <assert.h>
#include <unistd.h>

#ifndef PI
#define PI 3.14159265358979323846
#endif

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


typedef struct {
  int Np; // number of particles per string
  int Ns; // number of strings (cross sections) per cell
  double Hp; // length of string relative to cell size
  double eps; // threshold for convergence
  int itermax; // maximum number of iterations
  double eta; // relaxation factor
  bool nohf; // disable height fuctions, always use particles
  bool dumpcsv; // write csv files with particles and attraction points
  int stat_sh; // number of cells with height functions
  int stat_sc; // number of cells with particles
} Partstr;

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

const int kMaxSection = 125 * 2; // maximum number of endpoints
                                 // from neighbor cells
                                 // (5x5x5 stencil and 2 endpoints in each)
const int kMaxFacet = 12; // maximum number of vertices per facet

const int kMaxNp = 31; // maximum value of kPartstr.Np

static Trans kPartstrTrans; // if dumpcsv=1, contains last transformation
                            // used by GetCrossCurv()

#if dimension == 2
#define kNs 1
#else
#define kNs 3
#endif
static Partstr kPartstr = {7, kNs, 4., 1e-5, 20, 0.5, false, false};
#undef kNs

// Unit vector at angle ph.
static coord Unit(double ph) {
  coord p;
  p.x = cos(ph);
  p.y = sin(ph);
  p.z = 0;
  return p;
}

// Rotate planar vector e by planar unit vector de
static coord Rotate(coord e, coord de) {
  coord p;
  p.x = e.x * de.x - e.y * de.y;
  p.y = e.x * de.y + e.y * de.x;
  p.z = 0;
  return p;
}

// Rotate planar vector e by planar unit vector de in the negative direction
static coord Rotatem(coord e, coord de) {
  coord p;
  p.x = e.x * de.x + e.y * de.y;
  p.y = -e.x * de.y + e.y * de.x;
  p.z = 0;
  return p;
}

// Third component of cross product
static double Cross3(coord a, coord b) {
  return a.x * b.y - a.y * b.x;
}

static coord Cross(coord a, coord b) {
  coord r;
  r.x = a.y * b.z - a.z * b.y;
  r.y = a.z * b.x - a.x * b.z;
  r.z = a.x * b.y - a.y * b.x;
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

static double Norm1(coord a) {
  return fabs(a.x) + fabs(a.y) + fabs(a.z);
}

static double NormMax(coord a) {
  double r = fabs(a.x);
  if (fabs(a.y) > r) {
    r = fabs(a.y);
  }
  if (fabs(a.z) > r) {
    r = fabs(a.z);
  }
  return r;
}

static double Sqdist(coord a, coord b) {
  return Sqnorm(Sub(a, b));
}

static double Dist(coord a, coord b) {
  return sqrt(Sqdist(a, b));
}

static double Dotv(int np, const coord* aa, const coord* bb) {
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
  coord ep = Mul(Unit(ph + 0.5 * th), hp);
  coord em = Mul(Unit(ph - 0.5 * th), -hp);
  coord de = Unit(th);
  for (int j = 0; j < c; ++j) {
    xx[c + j + 1] = Add(xx[c + j], ep);
    xx[c - j - 1] = Add(xx[c - j], em);
    ep = Rotate(ep, de);
    em = Rotatem(em, de);
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
static void DxDph(coord p, double ph, double th, int np, double hp, coord* xx) {
  (void)p;
  int c = np / 2;
  xx[c] = Zero();
  coord ep = Mul(Unit(ph + 0.5 * th + PI * 0.5), hp);
  coord em = Mul(Unit(ph - 0.5 * th + PI * 0.5), -hp);
  coord de = Unit(th);
  for (int j = 0; j < c; ++j) {
    xx[c + j + 1] = Add(xx[c + j], ep);
    xx[c - j - 1] = Add(xx[c - j], em);
    ep = Rotate(ep, de);
    em = Rotatem(em, de);
  }
}

// Derivative dX/dth
static void DxDth(coord p, double ph, double th, int np, double hp, coord* xx) {
  (void)p;
  int c = np / 2;
  xx[c] = Zero();
  coord ep = Mul(Unit(ph + 0.5 * th + PI * 0.5), hp);
  coord em = Mul(Unit(ph - 0.5 * th + PI * 0.5), hp);
  coord de = Unit(th);
  for (int j = 0; j < c; ++j) {
    double jp = j + 0.5;
    xx[c + j + 1] = Add(xx[c + j], Mul(ep, jp));
    xx[c - j - 1] = Add(xx[c - j], Mul(em, jp));
    ep = Rotate(ep, de);
    em = Rotatem(em, de);
  }
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
static void F(
    int np, const coord* xx, int nl, const coord* ll, double eta, double k,
    coord* ff) {
  if (nl == 0) {
    for (int i = 0; i < np; ++i) {
      ff[i] = Zero();
    }
    return;
  }
  for (int i = 0; i < np; ++i) {
    coord pm = Nearest(ll[0], ll[1], xx[i]);

    for (int l = 0; l < nl; l += 2) {
      coord p = Nearest(ll[l], ll[l + 1], xx[i]);
      if (Sqdist(xx[i], p) < Sqdist(xx[i], pm)) {
        pm = p;
      }
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
static double Iter(
    coord* p_, double* ph_, double* th_, int np, double hp, coord* ff,
    coord* xx) {
  coord p = *p_;
  double ph = *ph_;
  double th = *th_;

  coord tt[kMaxNp];
  int c = np / 2;

  X(p, ph, th, np, hp, xx);

  // Skip correction of p, assume that p is initialized on a line segment
  // // correct p
  // p = Add(p, ff[c]);
  // X(p, ph, th, np, hp, tt);
  // for (int i = 0; i < np; ++i) {
  //   coord dx = Sub(tt[i], xx[i]);
  //   xx[i] = Add(xx[i], dx);
  //   ff[i] = Sub(ff[i], dx);
  // }

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

static double Curv(double hp, double th) {
  return 2. * sin(th * 0.5) / hp;
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

// Sets component i of point p to a
static void Set(coord* p, int i, double a) {
  if (i == 0) {
    p->x = a;
  } else if (i == 1) {
    p->y = a;
  } else {
    p->z = a;
  }
}

// Returns 3D base aligned with mesh directions.
// n: unit vector
// Output:
// *t, *u: vectors to form orthonormal base <t,u,n>
static void GetBase(coord n, coord* t, coord* u) {
  int i = Argmin(Abs(n));
  coord e = Zero();
  Set(&e, i, 1);
  n = Div(n, Norm(n));
  *t = Cross(n, e);
  *t = Div(*t, Norm(*t));
  *u = Cross(*t, n);
}

// Transformation from global to local coordinates.
// p: global coordinates
// o: origin
// t,n,u: orthonormal base
// Output:
// *l: local coordinates of p in <o,t,n,u>
static coord GlbToLoc(coord p, const Trans* w) {
  p = Sub(p, w->o);
  coord l = {Dot(w->t, p), Dot(w->n, p), Dot(w->u, p)};
  return l;
}

// Transformation from local to global coordinates.
// l: local coordinates in <t,n,u,o>
// t,n,u: orthonormal base
// o: origin
// Output:
// *p: global coordinates
static coord LocToGlb(coord l, const Trans* w) {
  coord p = w->o;
  p = Add(p, Mul(w->t, l.x));
  p = Add(p, Mul(w->n, l.y));
  p = Add(p, Mul(w->u, l.z));
  return p;
}

// Set z-copmonent to zero if 2D.
static void SetZeroZ(coord* p) {
  (void)p;
#if dimension == 2
  p->z = 0;
#endif
}

static int Facets(coord m, double alpha, coord* pp) {
#if dimension == 2
  int nf = facets(m, alpha, pp);
  for (int i = 0; i < nf; ++i) {
    SetZeroZ(&pp[i]);
  }
  return nf;
#else
  return facets(m, alpha, pp, 1.);
#endif
}

static coord Mycs(Point point, scalar c) {
  coord m = mycs(point, c);
  SetZeroZ(&m);
  return m;
}

// Appends csv file.
// incid: if true increment id
// c: attribute
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

// Returns the number of interfacial cells.
static int GetNcInter(scalar c) {
  int nc = 0;
  foreach () {
    if (interfacial(point, c)) {
      ++nc;
    }
  }
  return nc;
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
      double alpha = plane_alpha (c[], m);
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
  WriteVtkPoly(fn, xx, nx, pp, np, ss, "comment", true);

  free(ss);
  free(pp);
  free(xx);
}

// Cross-section of 2D interface from neighbor cells in plane coordinates.
// point: center cell
// c: volume fraction
// w: transformation
// ll: buffer for at least kMaxSection more points
// *nl: current size of ll
// Output:
// ll: appended by local coordinates of endpoints,
//     p = o + t*l.x  + t*l.y ,  [pa,pb] is one line segment
// *nl: new size of ll
static void Section2(
    Point point, scalar c, vector nn, const Trans* w, coord* ll, int* nl) {
  foreach_neighbor(2) {
    if (c[] > 0. && c[] < 1.) {
      coord m = {nn.x[], nn.y[], nn.z[]};
      double alpha = plane_alpha(c[], m);
      coord pp[kMaxFacet];
      int nf = Facets(m, alpha, pp);
      coord rn = {x, y, z};
      if (nf == 2) {
        coord p = Add(rn, Mul(pp[0], Delta));
        coord pb = Add(rn, Mul(pp[1], Delta));
        coord l = GlbToLoc(p, w);
        coord lb = GlbToLoc(pb, w);

        coord dl = Sub(l, lb);
        double mt = Dot(w->t, m);
        double mn = Dot(w->n, m);
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
static void Section3(
    Point point, scalar c, vector nn, const Trans* w, coord* ll, int* nl) {
  foreach_neighbor(2) {
    if (c[] > 0. && c[] < 1.) {
      int q = 0; // number of intersections found
      coord rn = {x, y, z}; // neighbor cell center

      // if cell intersects plane
      if (fabs(Dot(w->u, Sub(w->o, rn))) <= Delta * Norm1(w->u) * 0.5) {
        coord m = {nn.x[], nn.y[], nn.z[]}; // normal to facet
        double alpha = plane_alpha(c[], m); // plane constant
        coord pp[kMaxFacet];
        int nf = Facets(m, alpha, pp);
        assert(nf <= kMaxFacet);

        // find two intersections with facet edge and plane
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
              double mt = Dot(w->t, m); // local coordinates of normal
              double mn = Dot(w->n, m);
              double dt = dl.x; // local coordinates of line segment
              double dn = dl.y;
              // ensure positive orientation with normal m
              if (Cross3(Coord(dt, dn, 0), Coord(mt, mn, 0)) < 0) {
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
  }
}

// Cross-section of interface from neighbor cells in plane coordinates.
// point: center cell
// c: volume fraction
// w: transformation
// ll: buffer for at least kMaxSection more points
// *nl: current size of ll
// Output:
// ll: appended by local coordinates of endpoints,
//     p = o + t*l.x  + t*l.y ,  [pa,pb] is one line segment
// *nl: new size of ll
static void Section(
    Point point, scalar c, vector nn, const Trans* w, coord* ll, int* nl) {
#if dimension == 2
  Section2(point, c, nn, w, ll, nl);
#else
  Section3(point, c, nn, w, ll, nl);
#endif
}

// Curvature of a set line segments.
// ll: flat array endpoints of line segments
// nl: size of ll
// delta: cell size
// Output:
// res_: difference at last iteration
// it_: number of iterations
static double GetLinesCurv(
    coord* ll, int nl, double delta, const Partstr* conf, double* res_,
    int* it_) {
  if (nl >= 4) { // require at least two segments
    const int Np = conf->Np;
    const double eta = conf->eta;
    const int itermax = conf->itermax;
    const double hp = conf->Hp * delta / (Np - 1);

    coord xx[kMaxNp]; // positions
    coord ff[kMaxNp]; // forces
    coord p = {0., 0, 0.};
    double ph = 0.;
    double th = 0.;
    double k = 0;
    double res = 0;
    int it = 0;
    for (it = 0; it < itermax; ++it) {
      X(p, ph, th, Np, hp, xx);
      F(Np, xx, nl, ll, eta, k, ff);

      k = Curv(hp, th);
      res = Iter(&p, &ph, &th, Np, hp, ff, xx);
      if (res / (eta * delta) < conf->eps) {
        break;
      }
    }

    *res_ = res;
    *it_ = it;
    if (conf->dumpcsv) {
      coord tt[kMaxNp];
      for (int i = 0; i < Np; ++i) {
        tt[i] = Add(xx[i], Mul(ff[i], 1. / eta));
        tt[i] = LocToGlb(tt[i], &kPartstrTrans);
      }
      for (int i = 0; i < Np; ++i) {
        xx[i] = LocToGlb(xx[i], &kPartstrTrans);
      }

      AppendCsv(t * 100, Np, xx, "part", 0, -k);
      AppendCsv(t * 100, Np, tt, "attr", 0, -k);
    }
    return k;
  }
  return 0;
}

// Curvature of cross section by plane through a.n,a.t and point a.o
// point: target cell
// c: volume fraction
// nn: normal
// w: transformation to local coordinates
static double GetCrossCurv(
    Point point, scalar c, vector nn, const Trans* w, const Partstr* conf) {
  coord ll[kMaxSection];
  int nl = 0;
  Section(point, c, nn, w, ll, &nl);
  double res;
  int it;
  if (conf->dumpcsv) {
    kPartstrTrans = *w;
  }
  return GetLinesCurv(ll, nl, Delta, conf, &res, &it);
}

// Transformation b rotated at angle.
// s: index of cross section
// Ns: number of cross sections
static Trans GetSectionTrans(int s, int Ns, const Trans* b) {
  const double g = PI * s / Ns;
  Trans w = *b;
  w.t = Add(Mul(b->t, cos(g)), Mul(b->u, sin(g)));
  w.u = Cross(w.t, w.n);
  return w;
}

// Mean curvature over multiple cross sections by planes rotated around b.n
static double GetMeanCurv(
    Point point, scalar c, vector nn, const Trans* b, const Partstr* conf) {
  double ksum = 0;
  const int Ns = conf->Ns;
  for (int s = 0; s < Ns; ++s) {
    const Trans w = GetSectionTrans(s, Ns, b);
    ksum += GetCrossCurv(point, c, nn, &w, conf);
  }
  return ksum / Ns;
}

// Returns transformation to local coordinates at the interface.
static Trans GetPointTrans(Point point, scalar c, vector nn) {
  Trans b;
  coord n = {nn.x[], nn.y[], nn.z[]};
  double alpha = plane_alpha(c[], n);

  coord o;
  plane_area_center(n, alpha, &o);
  SetZeroZ(&o);
  coord rc = {x, y, z};
  o = Add(rc, Mul(o, Delta));

  b.o = o;
  b.n = Div(n, Norm(n));
  GetBase(b.n, &b.t, &b.u);
  return b;
}

static void CalcNormal(scalar c, vector nn) {
  foreach () {
    coord n = Mycs(point, c);
    nn.x[] = n.x;
    nn.y[] = n.y;
    nn.z[] = n.z;
  }
}

// Dumps cross sections of the inteface fragments to vtk.
void DumpLines(scalar c, vector nn, const Partstr* conf, const char* fn) {
  const int Ns = conf->Ns;

  const int nc = GetNcInter(c); // number of interfacial cells
  const int mm = nc * kMaxSection;

  coord* xx = (coord*)malloc(mm * sizeof(coord));
  int* pp = (int*)malloc(mm * sizeof(int));
  int* ss = (int*)malloc(mm * sizeof(int));
  int nx = 0;
  int np = 0;

  foreach() {
    if (interfacial(point, c)) {
      const Trans b = GetPointTrans(point, c, nn);

      for (int s = 0; s < Ns; ++ s) {
        const Trans w = GetSectionTrans(s, Ns, &b);

        coord ll[kMaxSection];
        int nl = 0;

        Section(point, c, nn, &w, ll, &nl);
        GetCrossCurv(point, c, nn, &w, conf);

        for (int i = 0; i < nl; ++i) {
          ll[i] = LocToGlb(ll[i], &w);
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
  WriteVtkPoly(fn, xx, nx, pp, np, ss, "lines", false);

  free(ss);
  free(pp);
  free(xx);
}


static double partstr(Point point, scalar c, vector nn, const Partstr* conf) {
  const Trans b = GetPointTrans(point, c, nn);
  double k = -GetMeanCurv(point, c, nn, &b, conf);
#if dimension == 3
  k *= 2;
#endif
  return k;
}

void CheckConf(const Partstr* conf) {
  if (!(conf->Np <= kMaxNp)) {
    fprintf(
        stderr, "error: Too many particles per string Np=%d > %d\n", conf->Np,
        kMaxNp);
    exit(1);
  }
  if (!(conf->Np % 2 == 1)) {
    fprintf(stderr, "error: Np=%d is not an odd number\n", conf->Np);
    exit(1);
  }
}

trace cstats curvature_partstr(struct Curvature p) {
  scalar c = p.c, kappa = p.kappa;
  double sigma = p.sigma ? p.sigma : 1.;
  int sh = 0, sf = 0, sa = 0, sc = 0;
  vector ch = c.height, h = automatic(ch);
  if (!ch.x.i) heights(c, h);

  Partstr* conf = &kPartstr;
  CheckConf(conf);

#if TREE
  kappa.refine = kappa.prolongation = curvature_prolongation;
  kappa.restriction = curvature_restriction;
#endif

  scalar k[];
  scalar_clone(k, kappa);

  vector nn[];
  CalcNormal(c, nn);

  boundary({nn});

  foreach (reduction(+:sh) reduction(+:sc)) {
    if (!interfacial(point, c)) {
      k[] = nodata;
    } else if (
        !conf->nohf && (k[] = height_curvature(point, c, h)) != nodata) {
      sh++;
    } else {
      k[] = partstr(point, c, nn, conf);
      sc++;
    }
  }
  boundary({k});

  foreach () {
    double kf = k[];
    if (kf == nodata)
      kappa[] = nodata;
    else if (p.add)
      kappa[] += sigma * kf;
    else
      kappa[] = sigma * kf;
  }
  boundary({kappa});

  conf->stat_sh = sh;
  conf->stat_sc = sc;
  return (cstats){sh, sf, sa, sc};
}

trace cstats curvature_orig(struct Curvature p) {
  return curvature(p);
}
#define curvature curvature_partstr
