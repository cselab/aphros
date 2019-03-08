#undef NDEBUG
#include "geom/mesh.h"

// Returns true if a < b (lex starting from end)
template <class T, size_t d>
bool Cmp(const GVect<T, d>& a, const GVect<T, d>& b) {
  int i = a.size();
  while (i--) {
    if (a[i] != b[i]) {
      return a[i] < b[i];
    }
  }
  return false;
}

const int dim = 3;
using MIdx = GMIdx<dim>;
using Dir = GDir<dim>;
using Scal = double;
using Vect = GVect<Scal, dim>;

void TestMesh() {
  Rect<Vect> dom(Vect(0), Vect(1));
  using M = MeshStructured<Scal, dim>;
  MIdx b(0);
  MIdx s(5, 5, 1);
  int hl = 2;
  Vect doms = dom.GetDimensions();
  Vect h = dom.GetDimensions() / Vect(s);
  M m = InitUniformMesh<M>(dom, b, s, hl, true, s);

  // Total volume
  Scal v = 0.;
  for (auto i : m.Cells()) {
    v += m.GetVolume(i);
  }
  // Comm
  FieldCell<Scal> fc;
}

int main() {
  TestMesh();
}
