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
  MIdx s(2, 2, 1);    // size in cells
  int hl = 0;         // halos 
  return InitUniformMesh<M>(dom, b, s, hl, true, s);
}

void TestAvg() {
  PF;

  M m = GetMesh(); 
  solver::Sphavg<M> sa(m, 3);
}

int main() {
  TestAvg();
}
