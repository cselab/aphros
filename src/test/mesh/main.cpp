// Created by Petr Karnakov on 30.05.2018
// Copyright 2018 ETH Zurich

#undef NDEBUG
#include <array>
#include <cassert>
#include <iostream>
#include <sstream>

#include "geom/block.h"
#include "geom/mesh.h"
#include "geom/mesh.ipp"
#include "util/format.h"
#include "util/suspender.h"

// Returns true if a < b (lex starting from end)
template <class T, size_t d>
bool Less(const generic::Vect<T, d>& a, const generic::Vect<T, d>& b) {
  int i = a.size();
  while (i--) {
    if (a[i] != b[i]) {
      return a[i] < b[i];
    }
  }
  return false;
}

const int dim = DIM;
using MIdx = generic::MIdx<dim>;
using Dir = GDir<dim>;
using Scal = double;
using Vect = generic::Vect<Scal, dim>;

void TestBlock() {
  std::cout << __func__ << std::endl;
  const size_t halos = 1;
  MIdx oi(0); // origin inner
  MIdx si(2); // size inner
  MIdx oa = oi - MIdx(halos); // origin all
  MIdx sa = si + MIdx(2 * halos); // size all

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
    assert(Less(xp, x));
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
  std::cout << __func__ << std::endl;
  const Rect<Vect> dom(Vect(0), Vect(1));
  using M = MeshCartesian<Scal, dim>;
  const MIdx begin(0);
  const MIdx size(2);
  const int halos = 2;
  Vect doms = dom.GetDimensions();
  const Vect h = dom.GetDimensions() / Vect(size);
  M m{begin, size, dom, halos, true, true, size, 0};

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
  auto sh = size + MIdx(halos * 2); // size with halos
  PCMP(m.GetAllBlockCells().size(), sh.prod());
  size_t nf = 0;
  for (int q = 0; q < dim; ++q) {
    auto w = sh;
    ++w[q];
    nf += w.prod();
  }
  PCMP(m.GetAllBlockFaces().size(), nf);
  PCMP(m.GetAllBlockNodes().size(), (sh + MIdx(1)).prod());

  // Distance between centers
  for (auto i : m.Cells()) {
    Vect xi = m.GetCenter(i);
    for (auto q : m.Nci(i)) {
      Dir d(q.raw() / 2);
      Scal k = (q.raw() % 2 == 0 ? -1. : 1.);
      auto j = m.GetCell(i, q);
      Vect xj = m.GetCenter(j);
      CMP((xj - xi)[d.raw()], h[d.raw()] * k);
    }
  }

  // Index of opposite face
  for (auto d : m.dirs) {
    PCMP(m.GetOpposite(IdxNci(2 * d)).raw(), 2 * d + 1);
    PCMP(m.GetOpposite(IdxNci(2 * d + 1)).raw(), 2 * d);
  }

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
      s << index.GetMIdx(f) << index.GetDir(f).GetLetter();
      return s.str();
    };
    for (auto d : m.dirs) {
      const MIdx w = m.GetInBlockCells().GetBegin();
      const IdxFace f = index.GetIdx(w, Dir(d));
      std::cout << "f=" << str(f) << ":";
      for (auto q : m.Nci(IdxCell(0))) {
        std::cout << " " << str(m.GetFace(f, q));
      }
      std::cout << std::endl;
    }
  }

  {
    std::cout << "\nFace neighbor nodes\n";
    const auto& indexn = m.GetIndexNodes();
    const auto& indexf = m.GetIndexFaces();
    auto str = [&](IdxFace f) {
      std::stringstream s;
      s << indexf.GetMIdx(f) << indexf.GetDir(f).GetLetter();
      return s.str();
    };
    for (auto d : m.dirs) {
      const MIdx w = m.GetInBlockCells().GetBegin();
      const IdxFace f = indexf.GetIdx(w, Dir(d));
      std::cout << "f=" << str(f) << ":";
      for (size_t q = 0; q < m.kFaceNumNeighborNodes; ++q) {
        std::cout << " " << indexn.GetMIdx(m.GetNode(f, q));
      }
      std::cout << std::endl;
    }
  }

  {
    std::cout << "\nCell neighbor nodes\n";
    const auto& indexc = m.GetIndexCells();
    const auto& indexn = m.GetIndexNodes();
    const MIdx w = m.GetInBlockCells().GetBegin();
    const IdxCell c = indexc.GetIdx(w);
    std::cout << "c=" << indexc.GetMIdx(c) << ":";
    for (size_t q = 0; q < m.kCellNumNeighborNodes; ++q) {
      std::cout << " " << indexn.GetMIdx(m.GetNode(c, q));
    }
    std::cout << std::endl;
  }

  {
    std::cout << "\nCell neighbor cells\n";
    IdxCell c(0);
    for (auto q : m.Nci(c)) {
      IdxFace f = m.GetFace(c, q);
      PCMP(m.GetNci(c, f).raw(), q.raw());
    }
  }

  // Comm
  FieldCell<Scal> fc(m);
  m.Comm(&fc);
  auto p = dynamic_cast<typename M::CommRequestScal*>(m.GetComm()[0].get());
  CMP(p->field, &fc);

  // Stencil indices
  {
    std::cout << "\nStencilIndices" << std::endl;
    auto& indexc = m.GetIndexCells();
    IdxCell c(indexc.GetIdx(m.GetInBlockCells().GetBegin() + MIdx(1)));
    std::cout << "c: " << c.GetRaw() << " " << indexc.GetMIdx(c) << std::endl;
    for (IdxCell cn : m.Stencil(c)) {
      std::cout << cn.GetRaw() << " " << indexc.GetMIdx(cn) << std::endl;
    }
  }
}

void TestMeshIndices() {
  std::cout << __func__ << std::endl;
  const Rect<Vect> dom(Vect(0), Vect(1));
  using M = MeshCartesian<Scal, dim>;
  const MIdx begin(0);
  const MIdx size(2);
  const int halos = 0;
  M m{begin, size, dom, halos, true, true, size, 0};

  {
    std::cout << "\nm.GetAllBlockCells()" << std::endl;
    auto blockca = m.GetAllBlockCells();
    auto indexc = m.GetIndexCells();
    for (auto w : blockca) {
      std::cout << w << " " << indexc.GetIdx(w).GetRaw() << std::endl;
    }
    std::cout << std::endl;
  }

  {
    std::cout << "\nm.GetAllBlockFaces()" << std::endl;
    auto blockfa = m.GetAllBlockFaces();
    auto indexf = m.GetIndexFaces();
    for (auto p : blockfa) {
      std::cout << p.first << " " << p.second.GetLetter() << ","
                << indexf.GetIdx(p).GetRaw() << std::endl;
    }
    std::cout << std::endl;
  }

  {
    std::cout << "\n(MIdx,Dir) <-> IdxFace" << std::endl;
    GBlock<IdxFace, dim> blockf(begin, size);
    GIndex<IdxFace, dim> indexf(begin, size + MIdx(1));
    for (auto p : blockf) {
      auto f = indexf.GetIdx(p);
      auto pp = indexf.GetMIdxDir(f);
      auto ff = indexf.GetIdx(pp);
      std::cout << p.first << "," << p.second.GetLetter() << " " << f.GetRaw()
                << " | " << pp.first << "," << pp.second.GetLetter() << " "
                << ff.GetRaw() << std::endl;
      assert(f == ff);
      assert(p == pp);
    }
    std::cout << std::endl;

    using It = typename GBlock<IdxFace, dim>::iterator;
    It it(&blockf, MIdx(0), Dir(0));
    int count = 0;
    for (it = blockf.begin(); it != blockf.end(); ++it) {
      std::cout << (*it).first << " " << (*it).second.GetLetter() << std::endl;
      ++count;
    }
    assert(count == dim * size.prod() + (MIdx(size.prod()) / size).sum());
  }
}

#if DIM == 3
void TestNotation() {
  std::cout << __func__ << std::endl;
  const Rect<Vect> dom(Vect(0), Vect(1));
  using M = MeshCartesian<Scal, dim>;
  const MIdx begin(0);
  const MIdx size(2);
  const size_t halos = 1;
  const M m{begin, size, dom, halos, true, true, size, 0};

  {
    std::cout << "\nTestNotation: IdxCellMesh\n";
    for (auto c : m.CellsM()) {
      std::cout << NAMEVALUE(c) << " ";
      auto dx = m.direction(0);
      for (auto d : m.dirs) {
        (void)d;
        dx = dx.next();
        std::cout << c + dx << " ";
        std::cout << c - dx << " ";
      }
      std::cout << std::endl;
    }
    for (auto c : m.SuCellsM()) {
      std::cout << MIdx(c).sum() << ' ';
    }
    std::cout << std::endl;
    auto c = m(IdxCell());
    std::stringstream("0 0 1") >> c;
    std::cout << NAMEVALUE(c) << ' ' << NAMEVALUE(size_t(IdxCell(c))) << ' '
              << NAMEVALUE(m.GetCenter(c)) << ' ' << ' ' << NAMEVALUE(MIdx(c));
    std::cout << std::endl;
  }

  {
    std::cout << "\nTestNotation: IdxFaceMesh\n";
    for (auto f : m.FacesM()) {
      std::cout << NAMEVALUE(f) << " ";
      auto dx = m.direction(0);
      for (auto d : m.dirs) {
        (void)d;
        dx = dx.next();
        std::cout << f + dx << ' ';
        std::cout << f - dx << ' ';
      }
      std::cout << std::endl;
    }
    for (auto f : m.SuFacesM()) {
      std::cout << MIdx(f).sum() << ' ';
    }
    std::cout << std::endl;
    auto f = m(IdxFace());
    std::stringstream("0 1 1 y") >> f;
    auto c = m(IdxCell());
    std::stringstream("0 0 1") >> c;
    std::cout << NAMEVALUE(f) << '\n';
    std::cout << NAMEVALUE(size_t(IdxFace(f))) << '\n';
    std::cout << NAMEVALUE(f.direction()) << '\n';
    std::cout << NAMEVALUE(f.area) << '\n';
    std::cout << NAMEVALUE(f.surface()) << '\n';
    std::cout << NAMEVALUE(f.cell(0)) << ' ' << NAMEVALUE(f.cell(1)) << '\n';
    std::cout << NAMEVALUE(m.GetCenter(f)) << '\n';
    std::cout << NAMEVALUE(c) << '\n';
    std::cout << NAMEVALUE(f.cm()) << ' ' << NAMEVALUE(f.cp());
    std::cout << '\n';
    std::cout << NAMEVALUE(c.center()) << ' ' << c.center[0] << '\n';
    std::cout << NAMEVALUE(c.volume) << '\n';
    std::cout << NAMEVALUE(c.outward_factor(IdxNci(0))) << '\n';
    std::cout << NAMEVALUE(c.outward_surface(IdxNci(0))) << '\n';
    std::cout << NAMEVALUE(c.face(IdxNci(0))) << '\n';
    std::cout << NAMEVALUE(c.face(-m.direction(0))) << '\n';
    std::cout << NAMEVALUE(f.center()) << ' ' << f.center[0] << '\n';
    auto center = c.center();
    std::cout << NAMEVALUE(center) << ' ' << center[0] << '\n';
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
#else
void TestNotation() {}
#endif

int main() {
  TestBlock();
  TestMesh();
  TestNotation();
  TestMeshIndices();
}
