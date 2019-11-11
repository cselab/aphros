#undef NDEBUG
#include <sstream>
#include <iostream>
#include <cassert>

#include "geom/mesh.h"

const int dim = 3;
using MIdx = GMIdx<dim>;
using Dir = GDir<dim>;
using Scal = double;
using Vect = GVect<Scal, dim>;

#define CMP(a, b) \
  assert(Cmp(a, b)); 

// Print CMP
#define PCMP(a, b) \
  std::cerr << #a << "=" << a << ", " << #b << "=" << b << std::endl; \
  CMP(a, b); 

void TestCellColumn() {
  Rect<Vect> dom{Vect(0), Vect(1)};
  using M = MeshStructured<Scal, dim>;
  MIdx b(0); // lower index
  MIdx s(2);    // size in cells
  int hl = 2;         // halos
  Vect h = dom.GetDimensions() / Vect(s);
  M m = InitUniformMesh<M>(dom, b, s, hl, true, true, s, 0);
}

int main() {
  TestCellColumn();
}
