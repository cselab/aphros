#undef NDEBUG
#include <sstream>
#include <iostream>
#include <cassert>
#include <functional>
#include <cmath>

#include "geom/mesh.h"
#include "solver/solver.h"
#include "func/init_u.h"
#include "parse/vars.h"

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
  FieldCell<Scal> fcu(m);

  Vars par;
  par.SetStr("string", "init_vf", "circlels");
  par.SetStr("int", "dim", "3");
  par.SetStr("vect", "circle_c", "0.5 0.5 0.5");
  par.SetStr("double", "circle_r", "0.3");

  auto fu = CreateInitU<Mesh>(par);
  fu(fcu, m);

  std::cout << GetNcInter(fcu);
}
