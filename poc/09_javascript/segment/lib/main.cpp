#include <solver/normal.h>
#include <solver/reconst.h>
#include <geom/mesh.h>

#include "segment.h"

const int D = 5, dim = 3, edim = 2;
using MIdx = GMIdx<dim>;
using Dir = GDir<dim>;
using Vect = GVect<Scal, dim>;
using M = MeshStructured<Scal, dim>;
using R = Reconst<Scal>;
using U = solver::UNormal<M>;

struct Segment {
    Scal n[2*D*D], a[D*D], s[4*D*D];
} segment;

inline void Clip(Scal& a, Scal l, Scal u) {
  a = std::max(l, std::min(u, a));
}

inline void Clip(Scal& a) {
  Clip(a, 0., 1.);
}

static Scal GetLineA1(const GVect<Scal, 3>& n, Scal u) {
  Scal nx = std::abs(n[0]);
  Scal ny = std::abs(n[1]);
  if (ny < nx) {
    std::swap(nx, ny);
  }
  Clip(u);
  if (u < 0.5) {
      return R::GetLineA0(nx, ny, u);
  } else {
      return -R::GetLineA0(nx, ny, 1. - u);
  }
}

int segment_get(const Scal alpha[D*D], /**/ Scal **pn, Scal **pa, Scal **ps) {
    enum {X, Y, Z};
    Rect<Vect> dom(Vect(0), Vect(D));
    MIdx b(0);
    MIdx size(D, D, 1);
    int hl;
    FieldCell<Vect> fcn;
    FieldCell<Scal> fck;
    Vect u;
    MIdx w;
    Scal *n, *a, al, *s;
    int i;

    hl = 2;
    M m = InitUniformMesh<M>(dom, b, size, hl, true, size);
    FieldCell<Scal> fcu(m, 0);
    FieldCell<bool> fci(m, true);

    n = segment.n;
    a = segment.a;
    s = segment.s;

    i = 0;
    for (auto c : m.Cells()) {
      fcu[c] = alpha[i++];
    }

    U::CalcNormalYoung(m, fcu, fci, /**/ fcn);

    i = 0;
    for (auto c : m.Cells()) {
      u = fcn[c];
      n[i++] = u[X];
      n[i++] = u[Y];
    }

    i = 0;
    for (auto c : m.Cells()) {
        u = fcn[c];
        al = alpha[i];
        a[i] = GetLineA1(u, al);
        i++;
    }

    i = 0;
    for (auto c : m.Cells()) {
        u = fcn[c];
        al = alpha[i];
        auto seg = R::GetLineEnds(u, al, Vect(1));
        auto x = m.GetCenter(c);
        seg[0] += x;
        seg[1] += x;
        s[i++] = seg[0][X];
        s[i++] = seg[0][Y];
        s[i++] = seg[1][X];
        s[i++] = seg[1][Y];
    }

    *pn = n; *pa = a; *ps = s;
    return 0;
}

int segment_norm(int i, int j, const double *a, /**/ double *px, double *py) {
#define b(i, j) (a[(D)*(j) + (i)])
    double nx, ny, n;
    nx = (b(i+1,j+1)+2*b(i+1,j)+b(i+1,j-1)-b(i-1,j+1)-2*b(i-1,j)-b(i-1,j-1))/8;
    ny = (b(i+1,j+1)-b(i+1,j-1)+2*b(i,j+1)-2*b(i,j-1)+b(i-1,j+1)-b(i-1,j-1))/8;
    n =  -(fabs(nx) + fabs(ny));
    nx /= n;
    ny /= n;
    *px = nx; *py = ny;
    return 0;
#undef b
}
