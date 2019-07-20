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
#include "solver/partstrmesh.h"

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

Mesh _mesh = GetMesh(MIdx(16));

#include "chpartstr.h"

int main() {
  auto& m = _mesh;
  FieldCell<Scal> fcu(m);

  Vars par;
  /*
  par.SetStr("string", "init_vf", "circlels");
  par.SetStr("int", "dim", "3");
  par.SetStr("vect", "circle_c", "0.52 0.51 0.5");
  par.SetStr("double", "circle_r", "0.3");
  */

  par.SetStr("string", "init_vf", "list");
  par.SetStr("string", "list_path", "b.dat");
  par.SetStr("int", "list_ls", "2");
  par.SetStr("int", "dim", "3");

  auto fu = CreateInitU<Mesh>(par);
  fu(fcu, m);

  std::cout << GetNcInter(fcu);

  FieldCell<Scal> nx(m), ny(m), nz(m);
  vector nn = {nx, ny, nz};
  CalcNormal(fcu, nn);

  DumpFacets(fcu, "o.vtk");
  DumpLines(fcu, nn, kPartstr, "line.vtk");

  kPartstr.csv = true;

  for (auto c : m.Cells()) {
    if (fcu[c] > 0 && fcu[c] < 1) {
      std::cout << partstr(c, fcu, nn) << std::endl;
    }
  }

  auto& kp = kPartstr;
  auto ps = std::make_shared<typename PartStr<Scal>::Par>();
  ps->npmax = kp.Np;
  ps->relax = kp.eta;
  ps->leq = kp.Hp;

  auto pm = std::make_shared<typename solver::PartStrMesh<Mesh>::Par>();
  pm->tol = kp.eps;
  pm->ns = kp.Ns;
  pm->itermax = kp.itermax;
  pm->ps = ps;
  solver::PartStrMesh<Mesh> psm(m, pm);
}
