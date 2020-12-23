// Created by Petr Karnakov on 30.05.2018
// Copyright 2018 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <iostream>
#include <sstream>

#include "geom/mesh.h"
#include "geom/notation.h"

// Returns true if a < b (lex starting from end)
template <class T, size_t d>
bool Cmp(const generic::Vect<T, d>& a, const generic::Vect<T, d>& b) {
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
using Vect = generic::Vect<Scal, dim>;

void TestBlock() {
  const size_t hl = 1;
  MIdx oi(0); // origin inner
  MIdx si(2); // size inner
  MIdx oa = oi - MIdx(hl); // origin all
  MIdx sa = si + MIdx(2 * hl); // size all

  GBlockFaces<dim> bi(oi, si);
  GIndex<IdxFace, dim> ba(oa, sa);

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

#define CMP(a, b) assert(Cmp(a, b));

// Print CMP
#define PCMP(a, b)                                                    \
  std::cout << #a << "=" << a << ", " << #b << "=" << b << std::endl; \
  CMP(a, b);

void TestMesh() {
  Rect<Vect> dom(Vect(0., 1.5, 2.7), Vect(5.3, 4.1, 3.));
  using M = MeshStructured<Scal, dim>;
  MIdx begin(-2, -3, -4); // lower index
  MIdx size(5, 4, 3); // size in cells
  int hl = 2; // halos
  Vect doms = dom.GetDimensions();
  Vect h = dom.GetDimensions() / Vect(size);
  M m = InitUniformMesh<M>(dom, begin, size, hl, true, true, size, 0);

  // Total volume
  Scal v = 0.;
  for (auto i : m.Cells()) {
    v += m.GetVolume(i);
  }
  PCMP(v, doms.prod());

  // Cell volume
  {
    IdxCell c(0);
    PCMP(m.GetVolume(c), h.prod());
    for (auto i : m.AllCells()) {
      CMP(m.GetVolume(i), m.GetVolume(c));
    }
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
  auto sh = size + MIdx(hl * 2); // size with halos
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
      auto j = m.GetCell(i, n);
      Vect xj = m.GetCenter(j);
      CMP((xj - xi)[size_t(d)], h[size_t(d)] * k);
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

  {
    std::cout << "Cell neighbor cells\n";
    const auto& ic = m.GetIndexCells();
    for (auto c : m.Cells()) {
      std::cout << "c=" << ic.GetMIdx(c) << ":";
      for (auto q : m.Nci(c)) {
        std::cout << " " << ic.GetMIdx(m.GetCell(c, q));
      }
      std::cout << std::endl;
      break;
    }
  }

  {
    std::cout << "\nFace neighbor faces\n";
    const auto& index = m.GetIndexFaces();
    auto str = [&](IdxFace f) {
      std::stringstream s;
      s << index.GetMIdx(f) << "," << index.GetDir(f).GetLetter();
      return s.str();
    };
    for (auto d : {0, 1, 2}) {
      const IdxFace f = index.GetIdx(m.GetInBlockCells().GetBegin(), Dir(d));
      std::cout << "f=" << str(f) << ":";
      for (auto q : {0, 1, 2, 3, 4, 5}) {
        const auto fn = m.GetFace(f, q);
        std::cout << " " << str(fn);
      }
      std::cout << std::endl;
    }
  }

  {
    std::cout << "\nCell neighbor cells\n";
    IdxCell c(0);
    for (auto q : m.Nci(c)) {
      IdxFace f = m.GetFace(c, q);
      PCMP(m.GetNci(c, f), q);
      PCMP(m.GetNci(IdxCell(2), f), -1);
    }
  }

  // Comm
  FieldCell<Scal> fc;
  m.Comm(&fc);
  auto p = dynamic_cast<typename M::CommRequestScal*>(m.GetComm()[0].get());
  CMP(p->field, &fc);

  // Stencil indices
  {
    std::cout << "\nStencilIndices" << std::endl;
    auto& bc = m.GetIndexCells();
    IdxCell c(bc.GetIdx(m.GetInBlockCells().GetBegin() + MIdx(1)));
    std::cout << "c: " << c.GetRaw() << " " << bc.GetMIdx(c) << std::endl;
    for (IdxCell cn : m.Stencil(c)) {
      std::cout << cn.GetRaw() << " " << bc.GetMIdx(cn) << std::endl;
    }
  }
}

void TestNotation() {
  const Rect<Vect> dom(Vect(0), Vect(1));
  using M = MeshStructured<Scal, dim>;
  const MIdx begin(0);
  const MIdx size(1, 2, 3);
  const size_t halos = 2;
  const M m = InitUniformMesh<M>(dom, begin, size, halos, true, true, size, 0);

  {
    std::cout << "\nTestNotation: IdxCellMesh\n";
    for (auto c : m.CellsM()) {
      auto dx = m.direction(0);
      auto dxb = m.direction(0, 0);
      auto dy = dx >> 1;
      auto dyb = -dy;
      auto dz = dx >> 2;
      fassert_equal(dy, m.direction(1));
      fassert_equal(dz, m.direction(2));
      auto dzb = dz.orient(Vect(-1));
      std::cout << NAMEVALUE(c) << " ";
      std::cout << c + dx << " ";
      std::cout << c + dy << " ";
      std::cout << c + dz << " ";
      std::cout << c + dxb << " ";
      std::cout << c + dyb << " ";
      std::cout << c + dzb << " ";
      std::cout << c - dzb << " ";
      std::cout << std::endl;
    }
    for (auto c : m.SuCellsM()) {
      std::cout << MIdx(c).sum() << ' ';
    }
    std::cout << std::endl;
    auto c = m.cell();
    std::stringstream("0 0 1") >> c;
    std::cout << NAMEVALUE(c) << ' ' << NAMEVALUE(size_t(IdxCell(c))) << ' '
              << NAMEVALUE(m.GetCenter(c)) << ' ' << ' ' << NAMEVALUE(MIdx(c));
    std::cout << std::endl;
  }

  {
    std::cout << "\nTestNotation: IdxFaceMesh\n";
    for (auto f : m.FacesM()) {
      auto dx = m.direction(0);
      auto dy = dx >> 1;
      auto dz = dx >> 2;
      auto dxb = dx.orient(Vect(-1));
      auto dyb = dy.orient(Vect(-1));
      auto dzb = dz.orient(Vect(-1));
      std::cout << NAMEVALUE(f) << " ";
      std::cout << f + dx << " ";
      std::cout << f + dy << " ";
      std::cout << f + dz << " ";
      std::cout << f + dxb << " ";
      std::cout << f + dyb << " ";
      std::cout << f + dzb << " ";
      std::cout << f - dzb << " ";
      std::cout << std::endl;
    }
    for (auto f : m.SuFacesM()) {
      std::cout << MIdx(f).sum() << ' ';
    }
    std::cout << std::endl;
    auto f = m.face();
    std::stringstream("0 1 1 y") >> f;
    auto c = m.cell();
    std::stringstream("0 0 1") >> c;
    std::cout << NAMEVALUE(f) << ' ' << NAMEVALUE(size_t(IdxFace(f))) << ' '
              << NAMEVALUE(f.direction()) << ' ' << NAMEVALUE(m.GetCenter(f))
              << ' ' << NAMEVALUE(c) << ' ' << NAMEVALUE((c + f)) << ' '
              << NAMEVALUE((c - f));
    std::cout << std::endl;
    std::cout << NAMEVALUE(f.cm()) << ' ' << NAMEVALUE(f.cp()) << ' '
              << NAMEVALUE(m(f.cm)) << ' ' << NAMEVALUE(m(f.cp));
    std::cout << std::endl;
    std::cout << NAMEVALUE(c.center[0]) << ' ' << NAMEVALUE(c.center()) << ' '
              << NAMEVALUE(f.center[0]) << ' ' << NAMEVALUE(f.center());
    std::cout << std::endl;
  }

  {
    std::cout << "\nTestNotation: IdxCellMesh gradient\n";
    auto dx = m.direction(0);
    auto dy = m.direction(1);
    auto dz = m.direction(2);
    FieldCell<Scal> u(m);
    for (auto c : m.CellsM()) {
      u[c] = c.center[0];
    }
    FieldCell<Scal> v(m);
    for (auto c : m.CellsM()) {
      v[c] = u[c + dx] + u[c - dx] + u[c + dy] + u[c - dy] + u[c + dz] +
             u[c - dz] - u[c] * 6;
    }
    FieldFace<Scal> uf(m);
    for (auto f : m.FacesM()) {
      uf[f] = u[f.cm] + u[f.cp];
    }
  }
}

int main() {
  TestBlock();
  TestMesh();
  TestNotation();

  {
    Rect<Vect> dom(Vect(0.), Vect(1.));
    using M = MeshStructured<Scal, dim>;
    MIdx begin(1, 1, 1); // lower index
    MIdx size(2, 2, 2); // size in cells
    int hl = 0; // halos
    M m = InitUniformMesh<M>(dom, begin, size, hl, true, true, size, 0);

    std::cout << "\nm.GetAllBlockCells()" << std::endl;
    auto bca = m.GetAllBlockCells();
    auto bc = m.GetIndexCells();
    for (auto w : bca) {
      std::cout << w << " " << bc.GetIdx(w).GetRaw() << std::endl;
    }
    std::cout << std::endl;

    std::cout << "\nm.GetAllBlockFaces()" << std::endl;
    auto bfa = m.GetAllBlockFaces();
    auto bf = m.GetIndexFaces();
    for (auto p : bfa) {
      std::cout << p.first << " " << p.second.GetLetter() << ","
                << bf.GetIdx(p).GetRaw() << std::endl;
    }
    std::cout << std::endl;

    std::cout << "\n(MIdx,Dir) <-> IdxFace" << std::endl;
    GBlock<IdxFace, 3> sb(begin, size);
    GIndex<IdxFace, 3> ind(begin, size + MIdx(1));
    for (auto p : sb) { // p: (w,d)
      auto j = ind.GetIdx(p);
      auto pp = ind.GetMIdxDir(j);
      auto jj = ind.GetIdx(pp);
      std::cout << p.first << "," << p.second.GetLetter() << " " << j.GetRaw()
                << " | " << pp.first << "," << pp.second.GetLetter() << " "
                << jj.GetRaw() << std::endl;
      assert(j == jj);
      assert(p == pp);
    }
    std::cout << std::endl;

    {
      using I = typename GBlock<IdxFace, 3>::iterator;
      I i(&sb, MIdx(0, 0, 0), Dir(0));
      auto l = [&]() {
        std::cout << (*i).first << " " << (*i).second.GetLetter() << std::endl;
      };

      int c = 0;
      for (i = sb.begin(); i != sb.end(); ++i) {
        l();
        ++c;
      }
      assert(
          c == 3 * size.prod() + size[0] * size[1] + size[1] * size[2] +
                   size[2] * size[0]);
    }
  }
}
