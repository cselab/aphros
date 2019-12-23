#undef NDEBUG
#include "geom/mesh.h"
#include "solver/normal.h"

const int dim = 3, edim = 2;
using MIdx = GMIdx<dim>;
using Dir = GDir<dim>;
using Scal = double;
using Vect = GVect<Scal, dim>;
using M = MeshStructured<Scal, dim>;
using std::cout;

int main() {
  enum { X, Y, Z };
  Rect<Vect> dom(Vect(0), Vect(1));

  MIdx b(0);
  MIdx s(5, 5, 1);
  int hl = 2;
  M m = InitUniformMesh<M>(dom, b, s, hl, true, true, s, 0);
  FieldCell<Scal> fcu(m, 0);
  FieldCell<bool> fci(m, true);
  FieldCell<Vect> fcn;
  FieldCell<Scal> fck;
  Vect x;
  MIdx w;

  for (auto c : m.Cells()) {
    x = m.GetCenter(c);
    fcu[c] = x[X] * x[Y];
  }
  solver::UNormal<M>::CalcNormal(m, fcu, fci, edim, fcn);

  auto& bc = m.GetIndexCells();
  for (auto c : m.Cells()) {
    w = bc.GetMIdx(c);
    if (w == s / 2) {
      cout << w << '\n';
      cout << fcn[c] << '\n';
    }
  }

  //  cout << fcu;
}
