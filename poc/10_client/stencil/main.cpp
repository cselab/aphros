#include <geom/mesh.h>

const int dim = 3;
using MIdx = GMIdx<dim>;
using Dir = GDir<dim>;
using Scal = double;
using Vect = GVect<Scal, dim>;
using M = MeshStructured<Scal, dim>;
using std::cout;

int main() {
  enum {X, Y, Z};
  Rect<Vect> dom(Vect(0), Vect(1));

  MIdx b(0); // mesh origin
  MIdx s(5); // mesh size
  int hl = 2;
  M m = InitUniformMesh<M>(dom, b, s, hl, true, true, s, 0);
  FieldCell<Scal> fcu(m, 0);
  Vect x;
  MIdx w;

  for (auto c : m.Cells()) {
      x = m.GetCenter(c);
      fcu[c] = x[X]*x[Y];
  }

  auto& bc = m.GetIndexCells();

  const int sw = 1;   // stencil halfwidth
  // block of offsets
  GBlock<IdxCell, dim> bo(MIdx(-sw), MIdx(sw * 2 + 1)); 

  MIdx wc = s / 2;
  std::cout << "neighbors of cell wc=" << wc << std::endl;
  for (auto c : m.Cells()) {
    MIdx w = bc.GetMIdx(c);
    if (w == wc) {
      for (MIdx wo : bo) {
        MIdx wn = w + wo;
        IdxCell cn = bc.GetIdx(wn);
        std::cout << "wn=" << wn << " u[wn]=" << fcu[cn] << std::endl;
      }
    }
  }
}
