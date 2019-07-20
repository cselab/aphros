#include "solver/reconst.h"
#include "solver/normal.h"

struct coord {
  double x, y, z;
};

using scalar = FieldCell<Scal>&;

struct vector {
  FieldCell<Scal>& x;
  FieldCell<Scal>& y;
  FieldCell<Scal>& z;
};

bool interfacial(IdxCell point, scalar c) {
  return c[point] > 0 && c[point] < 1;
}

#define static 
#define dimension 3
#define Point IdxCell
#define CELL point

#define POINTXYZ \
  double x = (*_mesh).GetCenter(point)[0]; \
  double y = (*_mesh).GetCenter(point)[1]; \
  double z = (*_mesh).GetCenter(point)[2]; \
  (void) x; (void) y; (void) z;


#define foreach() for (auto point : (*_mesh).Cells()) { \
  POINTXYZ

#define foreach_end }

auto _bc = (*_mesh).GetIndexCells();

#define foreach_neighbor(W) { \
  GBlock<IdxCell, dim> _bo(MIdx(-W), MIdx(W * 2 + 1)); \
  MIdx _w = _bc.GetMIdx(point); \
  for (MIdx _wo : _bo) { \
    IdxCell point = _bc.GetIdx(_w + _wo); \
    POINTXYZ

#define foreach_neighbor_end }}

using std::min;
using std::max;

double t;
double Delta = (*_mesh).GetCellSize()[0];

static double clamp(double a, double a0, double a1) {
  return a > a1 ? a1 : a < a0 ? a0 : a;
}

static double sq(double a) {
  return a * a;
}

double plane_alpha(double u, coord n) {
  return Reconst<Scal>::GetLineA(Vect(n.x, n.y, n.z), u, Vect(1));
}

void plane_area_center(coord n, double a, coord* o) {
  auto x = Reconst<Scal>::GetCenter(Vect(n.x, n.y, n.z), a, Vect(1));
  o->x = x[0];
  o->y = x[1];
  o->z = x[2];
}

FieldCell<bool> DetectInterface(const FieldCell<Scal>& fcu) {
  FieldCell<bool> fci((*_mesh), false);
  for (auto c : (*_mesh).AllCells()) {
    if (fcu[c] > 0. && fcu[c] < 1.) { fci[c] = true; }
  }
  return fci;
}

static coord mycs(Point point, scalar fcu) {
  static FieldCell<Vect> fcn;
  if (fcn.empty()) {
    solver::UNormal<Mesh>::CalcNormal(
        (*_mesh), fcu, DetectInterface(fcu), dim, fcn);
  }
  auto n = fcn[point];
  coord r = {n[0], n[1], n[2]};
  return r;
}

static int facets(coord n, double a, coord* pp, double h_=1) {
  (void) h_;
  std::vector<Vect> vv =
      Reconst<Scal>::GetCutPoly(Vect(0), Vect(n.x, n.y, n.z), a, Vect(1));
  int i = 0;
  for (auto& v : vv) {
    coord& p = pp[i];
    p.x = v[0];
    p.y = v[1];
    p.z = v[2];
    ++i;
  }
  return i;
}

#include "partstr.h"
