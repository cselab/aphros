#undef NDEBUG
#include "geom/mesh.h"

const int dim = 3;
using MIdx = GMIdx<dim>;
using Dir = GDir<dim>;
using Scal = double;
using Vect = GVect<Scal, dim>;

int main() {
  Rect<Vect> dom(Vect(0), Vect(1));
  using M = MeshStructured<Scal, dim>;
  MIdx b(0);
  MIdx s(5, 5, 1);
  int hl = 2;
  M m = InitUniformMesh<M>(dom, b, s, hl, true, s);

  // Total volume
  Scal v = 0.;
  for (auto i : m.Cells()) {
    v += m.GetVolume(i);
  }
  // Comm
  FieldCell<Scal> fc;
}
