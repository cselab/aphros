#undef NDEBUG
#include "solver/normal.h"
#include "geom/mesh.h"

const int dim = 3, edim = 2;
using MIdx = GMIdx<dim>;
using Dir = GDir<dim>;
using Scal = double;
using Vect = GVect<Scal, dim>;

int main() {
  enum {X, Y, Z};
  using std::cout;

  Rect<Vect> dom(Vect(0), Vect(1));
  using M = MeshStructured<Scal, dim>;
  MIdx b(0);
  MIdx s(5, 5, 1);
  int hl = 2;
  M m = InitUniformMesh<M>(dom, b, s, hl, true, s);
  FieldCell<Scal> fcu(m);
  FieldCell<bool> fci(m, true);
  FieldCell<Vect> fcn;
  FieldCell<Scal> fck;
  Vect x;

  for (auto c : m.AllCells()) {
      x = m.GetCenter(c);
      fcu[c] = x[X];
  }
  
  solver::UNormal<M>::CalcNormal(m, fcu, fci, edim, fcn, fck);
  cout << fcu;
}
