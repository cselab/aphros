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
  FieldCell<Scal> fcu(m, 0);
  FieldCell<bool> fci(m, true);
  FieldCell<Vect> fcn;
  FieldCell<Scal> fck;
  Vect x;
  MIdx w;

  for (auto c : m.Cells()) {
      x = m.GetCenter(c);
      fcu[c] = x[X]*x[Y];
  }
  solver::UNormal<M>::CalcNormal(m, fcu, fci, edim, fcn, fck);

  auto& bc = m.GetIndexCells();
  for (auto c : m.Cells()) {
      w = bc.GetMIdx(c);
      if (w == s/2) {
          cout << w << '\n';
          cout << fcn[c] << '\n';
      }
  }

  //  cout << fcu;


}
