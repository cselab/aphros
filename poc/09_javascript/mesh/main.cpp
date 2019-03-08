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

bool Cmp(Scal a, Scal b) {
  return std::abs(a - b) < 1e-12;
}

template <class T>
bool Cmp(T* a, T* b) {
  return a == b;
}

bool Cmp(size_t a, size_t b) {
  return a == b;
}

#define CMP(a, b) \
  assert(Cmp(a, b)); 

// Print CMP
#define PCMP(a, b) \
  std::cerr << #a << "=" << a << ", " << #b << "=" << b << std::endl; \
  CMP(a, b); 

void TestMesh() {
  Rect<Vect> dom(Vect(0., 1.5, 2.7), Vect(5.3, 4.1, 3.));
  using M = MeshStructured<Scal, dim>;
  MIdx b(-2, -3, -4); // lower index
  MIdx s(5, 4, 3);    // size in cells
  int hl = 2;         // halos 
  Vect doms = dom.GetDimensions();
  Vect h = dom.GetDimensions() / Vect(s);
  M m = InitUniformMesh<M>(dom, b, s, hl, true, s);

  // Total volume
  Scal v = 0.;
  for (auto i : m.Cells()) {
    v += m.GetVolume(i);
  }
  PCMP(v, doms.prod());

  // Cell volume
  IdxCell c(0);
  PCMP(m.GetVolume(c), h.prod());
  for (auto i : m.AllCells()) {
    CMP(m.GetVolume(i), m.GetVolume(c));
  }

  // Face area
  for (int q = 0; q < dim; ++q) {
    Dir d(q);
    IdxFace f(0);
    // Find any face with direction d
    for (IdxFace i : m.AllFaces()) {
      if (m.GetDir(i) == d) {
        f = i;
        break;
      }
    }
    assert(m.GetDir(f) == d);

    PCMP(m.GetArea(f), h.prod() / h[q]);
    for (auto i : m.AllFaces()) {
      if (m.GetDir(i) == d) {
        CMP(m.GetArea(i), m.GetArea(f));
      }
    }
  }

  // Number of elements
  auto sh = s + MIdx(hl * 2); // size with halos
  PCMP(m.GetAllBlockCells().size(), sh.prod());
  size_t nf = 0;
  for (int q = 0; q < 3; ++q) {
    auto w = sh;
    ++w[q];
    nf += w.prod();
  }
  PCMP(m.GetAllBlockFaces().size(), nf);
  PCMP(m.GetAllBlockNodes().size(), (sh + MIdx(1)).prod());

  // Distance between centers
  for (auto i : m.Cells()) {
    Vect xi = m.GetCenter(i);
    for (auto n : m.Nci(i)) {
      Dir d(n / 2); 
      Scal k = (n % 2 == 0 ? -1. : 1.);
      auto j = m.GetNeighbourCell(i, n);
      Vect xj = m.GetCenter(j);
      CMP((xj-xi)[size_t(d)], h[size_t(d)] * k);
    }
  }

  // Index of opposite face
  PCMP(m.GetOpposite(0), 1);
  PCMP(m.GetOpposite(1), 0);
  PCMP(m.GetOpposite(2), 3);
  PCMP(m.GetOpposite(3), 2);
  PCMP(m.GetOpposite(4), 5);
  PCMP(m.GetOpposite(5), 4);
  PCMP(m.GetOpposite(-1), -1);
  
  // Index of neighbour face
  {
    IdxCell c(0);
    for (auto q : m.Nci(c)) {
      IdxFace f = m.GetNeighbourFace(c, q);
      PCMP(m.GetNci(c, f), q);
      PCMP(m.GetNci(IdxCell(2), f), -1);
    }
  }

  // Comm
  FieldCell<Scal> fc;
  m.Comm(&fc);
  auto p = dynamic_cast<typename M::CoFcs*>(m.GetComm()[0].get());
  CMP(p->f, &fc);
}

int main() {
  TestMesh();

  {
    Rect<Vect> dom(Vect(0.), Vect(1.));
    using M = MeshStructured<Scal, dim>;
    MIdx b(1, 1, 1); // lower index
    MIdx s(2, 2, 2);    // size in cells
    int hl = 0;         // halos 
    M m = InitUniformMesh<M>(dom, b, s, hl, true, s);

    std::cout << "m.GetAllBlockCells()" << std::endl;
    auto bca = m.GetAllBlockCells();
    auto bc = m.GetIndexCells();
    for (auto w : bca) {
      std::cout << w << " " << bc.GetIdx(w).GetRaw() << std::endl;
    }
    std::cout << std::endl;

    std::cout << "m.GetAllBlockFaces()" << std::endl;
    auto bfa = m.GetAllBlockFaces();
    auto bf = m.GetIndexFaces();
    for (auto p : bfa) {
      std::cout 
          << p.first << " " << p.second.GetLetter() << ","
          << bf.GetIdx(p).GetRaw()
          << std::endl;
    }
    std::cout << std::endl;

    std::cout << "(MIdx,Dir) <-> IdxFace" << std::endl;
    GBlock<IdxFace, 3> sb(b, s);
    GIndex<IdxFace, 3> ind(b, s + MIdx(1));
    for (auto p : sb) { // p: (w,d)
      auto j = ind.GetIdx(p);
      auto pp = ind.GetMIdxDir(j);
      auto jj = ind.GetIdx(pp);
      std::cout 
          << p.first << "," << p.second.GetLetter() << " " << j.GetRaw() 
          << " | "
          << pp.first << "," << pp.second.GetLetter() << " " << jj.GetRaw()
          << std::endl;
      assert(j == jj);
      assert(p == pp);
    }
    std::cout << std::endl;

    {
      using I = typename GBlock<IdxFace, 3>::iterator;
      I i(&sb, MIdx(0, 0, 0), Dir::i);
      auto l = [&]() {
          std::cout << (*i).first << " " << (*i).second.GetLetter() << std::endl;
        };

      int c = 0;
      for (i = sb.begin(); i != sb.end(); ++i) {
        l();
        ++c;
      }
      assert(c == 3 * s.prod() + s[0] * s[1] + s[1] * s[2] + s[2] * s[0]); 
    }
  }
}
