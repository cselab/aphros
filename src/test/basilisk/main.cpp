#undef NDEBUG
#include <sstream>
#include <iostream>
#include <cassert>
#include <functional>
#include <cmath>

#include "geom/mesh.h"
#include "solver/solver.h"

const int dim = 3;
using MIdx = GMIdx<dim>;
using IdxCell = IdxCell;
using IdxFace = IdxFace;
using Dir = GDir<dim>;
using Scal = double;
using Vect = GVect<Scal, dim>;
using Mesh = MeshStructured<Scal, dim>;


Mesh GetMesh(MIdx s /*size in cells*/) {
  Rect<Vect> dom(Vect(0), Vect(1));
  MIdx b(0, 0, 0); // lower index
  int hl = 2;         // halos 
  return InitUniformMesh<Mesh>(dom, b, s, hl, true, true, s, 0);
}

Mesh _mesh = GetMesh(MIdx(8));

#include "chpartstr.h"

int main() {
  auto& m = _mesh;
  FieldCell<Scal> fc(m);
  for (auto c : m.Cells()) {
    fc[c] = (m.GetCenter(c).dist(Vect(0.5)) < 0.3 ? 0.5 : 0.);
  }
  std::cout << GetNcInter(fc);
}
