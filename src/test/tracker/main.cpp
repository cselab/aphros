#undef NDEBUG
#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <string>
#include <limits>
#include <stdexcept>

#include "debug/isnan.h"
#include "geom/field.h"
#include "geom/range.h"
#include "geom/mesh.h"
#include "solver/sphavg.h"

// try-catch
#define TR(...) try { __VA_ARGS__ } \
  catch (const std::runtime_error& e) { std::cout << e.what() << std::endl; }

// print current func
#define PF std::cout << std::endl << __func__ << std::endl;


const int dim = 3;
using MIdx = GMIdx<dim>;
using Dir = GDir<dim>;
using Scal = double;
using Vect = GVect<Scal, dim>;
using M = MeshStructured<Scal, dim>;

M GetMesh() {
  Rect<Vect> dom(Vect(0., 0., 0.), Vect(1., 1., 1.));
  using M = MeshStructured<Scal, dim>;
  MIdx b(0, 0, 0); // lower index
  MIdx s(10, 10, 10);    // size in cells
  int hl = 2;         // halos 
  return InitUniformMesh<M>(dom, b, s, hl, true, true, MIdx(20, 20, 20), 0);
}

void TestAvg() {
  PF;

  using Sphavg = solver::Sphavg<M>;

  M m = GetMesh(); 
  solver::Sphavg<M> sa(m, 2);

  FieldCell<Scal> fcu(m);
  FieldCell<Vect> fcv(m);
  FieldCell<Vect> fcvm(m);
  FieldCell<Scal> fcp(m);
  Scal dt = 1.;
  using Sph = typename Sphavg::Sph;
  std::vector<Sph> ss;
  ss.emplace_back(Vect(1.9), 0.2, 0.05);

  sa.Update(fcu, fcv, fcvm, dt, fcp, ss);
  sa.Update(fcu, fcv, fcvm, dt, fcp, ss);
  sa.Update(fcu, fcv, fcvm, dt, fcp, ss);
  sa.Update(fcu, fcv, fcvm, dt, fcp, ss);
}

int main() {
  TestAvg();
}
