#undef NDEBUG
#include <sstream>
#include <iostream>
#include <cassert>

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
using IdxCell = IdxCell;
using IdxFace = IdxFace;
using Dir = GDir<dim>;
using Scal = double;
using Vect = GVect<Scal, dim>;

void TestBlock() {
  const size_t hl = 1;
  MIdx oi(0); // origin inner
  MIdx si(2); // size inner
  MIdx oa= oi - MIdx(hl); // origin all
  MIdx sa = si + MIdx(2 * hl); // size all

  GBlockFaces<dim> bi(oi, si);
  GBlockFaces<dim> ba(oa, sa);

  GRange<IdxFace> ra(ba);
  GRangeIn<IdxFace, dim> ri(ba, bi);

  const MIdx xp0 = oa - MIdx(Dir(0));
  MIdx xp = xp0;
  Dir dp(0); // direction

  // Check that whole inner block covered with ascending indices
  for (auto i : ri) {
    auto x = ba.GetMIdx(i);
    auto d = ba.GetDir(i);

    // Next direction, reset xp
    if (dp < d) {
      xp = xp0;
    }

    // std::cerr << x << " " << d.GetLetter() << std::endl;
    assert(Cmp(xp, x));
    assert(oi <= x && x < oi + si + MIdx(d));

    xp = x;
    dp = d;
  }
}

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
  M m = InitUniformMesh<M>(dom, b, s, hl);

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
    IdxFace f;
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
  PCMP(m.GetNumCells(), sh.prod());
  size_t nf = 0;
  for (int q = 0; q < 3; ++q) {
    auto w = sh;
    ++w[q];
    nf += w.prod();
  }
  PCMP(m.GetNumFaces(), nf);
  PCMP(m.GetNumNodes(), (sh + MIdx(1)).prod());

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


  // Comm
  FieldCell<Scal> fc;
  m.Comm(&fc);
  auto p = dynamic_cast<typename M::CoFcs*>(m.GetComm()[0].get());
  CMP(p->f, &fc);
}

int main() {
  TestBlock();

  TestMesh();


  {
    Rect<Vect> dom(Vect(0.), Vect(1.));
    using M = MeshStructured<Scal, dim>;
    MIdx b(3, 2, 5); // lower index
    MIdx s(2, 2, 2);    // size in cells
    int hl = 0;         // halos 
    M m = InitUniformMesh<M>(dom, b, s, hl);
    for (auto i : m.GetBlockCells()) {
      std::cout << i << std::endl;
    }
    auto bf = m.GetBlockFaces();
    for (auto i : bf) {
      std::cout 
          << i.first << " " << i.second.GetLetter() << " "
          << bf.GetIdx(i).GetRaw() << " " 
          << std::endl;
    }
    std::cout << std::endl;

    using B = GBlock<IdxFace, 3>;
    B sb(b, s);
    for (auto i : sb) {
      auto j = sb.GetIdx(i);
      auto p = sb.GetMIdxDir(j);
      auto jj = sb.GetIdx(p);
      assert(j == jj);
      std::cout 
          << i.first << " " << i.second.GetLetter() << " | "
          << j.GetRaw() << " | "
          << p.first << " " << p.second.GetLetter() << " | "
          << jj.GetRaw() << " "
          << std::endl;
    }

    {
      using I = typename B::iterator;
      I i(&sb, MIdx(0, 0, 0), Dir::i);
      auto l = [&]() {
          std::cout << (*i).first << " " << (*i).second.GetLetter() << std::endl;
        };

      size_t c = 0;
      for (i = sb.begin(); i != sb.end(); ++i) {
        l();
        ++c;
      }
      assert(c == 3 * s.prod() + s[0] * s[1] + s[1] * s[2] + s[2] * s[0]); 
    }
  }

}
