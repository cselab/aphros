// Created by Petr Karnakov on 30.05.2018
// Copyright 2018 ETH Zurich

#include <cassert>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

#include "geom/mesh.h"
#include "solver/approx.h"
#include "solver/approx_eb.h"
#include "solver/cond.h"
#include "solver/solver.h"
#include "util/sysinfo.h"
#include "util/timer.h"

#define EXPOSE(x) do { \
  volatile auto EXPOSE_vol = x; \
  (void) EXPOSE_vol; \
} while (0);

const int dim = 3;
using MIdx = GMIdx<dim>;
using IdxCell = IdxCell;
using IdxFace = IdxFace;
using Dir = GDir<dim>;
using Scal = double;
using Vect = generic::Vect<Scal, dim>;
using M = MeshStructured<Scal, dim>;

// Echo Execute
#define EE(...)                                   \
  ;                                               \
  std::cerr << "\n" << #__VA_ARGS__ << std::endl; \
  __VA_ARGS__;

M GetMesh(MIdx size) {
  Rect<Vect> dom(Vect(0.1, 0.2, 0.1), Vect(1.1, 1.2, 1.3));
  MIdx begin(-2, -3, -4);
  int halos = 2;
  return InitUniformMesh<M>(dom, begin, size, halos, true, true, size, 0);
}

class TimerMesh : public Timer {
 public:
  TimerMesh(const std::string& name, M& m_) : Timer(name), m(m_) {}

 protected:
  M& m;
};

/*
class IFactoryTimerMesh {
 public:
  virtual TimerMesh* Make(M& m) = 0;
};

std::vector<std::unique_ptr<IFactoryTimerMesh>> reg;

template <class T>
class FactoryTimerMesh : IFactoryTimerMesh {
 public:
  TimerMesh* Make(M& m) override {
    return new T(m);
  }
};
*/

class Empty : public TimerMesh {
 public:
  Empty(M& m_) : TimerMesh("empty", m_) {}
  void F() override {}
};

class LoopPlain : public TimerMesh {
 public:
  LoopPlain(M& m_) : TimerMesh("loop-plain", m_) {}
  void F() override {
    size_t a = 0;
    for (size_t i = 0; i < m.GetAllBlockCells().size(); ++i) {
      a += i;
    }
    EXPOSE(a);
  }
};

class LoopAllCells : public TimerMesh {
 public:
  LoopAllCells(M& m_) : TimerMesh("loop-allcells", m_) {}
  void F() override {
    size_t a = 0;
    for (auto i : m.AllCells()) {
      a += i.GetRaw();
    }
    EXPOSE(a);
  }
};

class LoopInCells : public TimerMesh {
 public:
  LoopInCells(M& m_) : TimerMesh("loop-incells", m_) {}
  void F() override {
    size_t a = 0;
    for (auto i : m.Cells()) {
      a += i.GetRaw();
    }
    EXPOSE(a);
  }
};

class LoopAllFaces : public TimerMesh {
 public:
  LoopAllFaces(M& m_) : TimerMesh("loop-allfaces", m_) {}
  void F() override {
    size_t a = 0;
    for (auto i : m.AllFaces()) {
      a += i.GetRaw();
    }
    EXPOSE(a);
  }
};

class LoopInFaces : public TimerMesh {
 public:
  LoopInFaces(M& m_) : TimerMesh("loop-infaces", m_) {}
  void F() override {
    size_t a = 0;
    for (auto i : m.Faces()) {
      a += i.GetRaw();
    }
    EXPOSE(a);
  }
};

// Loop access to field
class LoopFldPlain : public TimerMesh {
 public:
  LoopFldPlain(M& m_)
      : TimerMesh("loop-fld-plain", m_), v(GRange<IdxCell>(m).size()) {}
  void F() override {
    Scal a = 0;
    for (size_t i = 0; i < v.size(); ++i) {
      v[i] += a;
      a = v[i];
    }
    EXPOSE(a);
  }

 private:
  std::vector<Scal> v;
};

class LoopFldAllCells : public TimerMesh {
 public:
  LoopFldAllCells(M& m_) : TimerMesh("loop-fld-allcells", m_), v(m) {}
  void F() override {
    Scal a = 0;
    for (auto i : m.AllCells()) {
      v[i] += a;
      a = v[i];
    }
    EXPOSE(a);
  }

 private:
  FieldCell<Scal> v;
};

class LoopMIdxAllCells : public TimerMesh {
 public:
  LoopMIdxAllCells(M& m_) : TimerMesh("loop-midx-allcells", m_) {}
  void F() override {
    size_t a = 0;
    for (auto c : m.AllCells()) {
      auto w = m.GetIndexCells().GetMIdx(c);
      a += w[0];
    }
    EXPOSE(a);
  }
};

class LoopMIdxAllFaces : public TimerMesh {
 public:
  LoopMIdxAllFaces(M& m_) : TimerMesh("loop-midx-allfaces", m_) {}
  void F() override {
    size_t a = 0;
    for (auto f : m.AllFaces()) {
      auto wd = m.GetIndexFaces().GetMIdxDir(f);
      a += wd.first[0];
      a += wd.second.GetLetter();
    }
    EXPOSE(a);
  }
};

class CellVolume : public TimerMesh {
 public:
  CellVolume(M& m_) : TimerMesh("m-cell-volume", m_) {}
  void F() override {
    size_t a = 1;
    for (auto c : m.Cells()) {
      a += m.GetVolume(c);
    }
    EXPOSE(a);
  }
};

class CellCenter : public TimerMesh {
 public:
  CellCenter(M& m_) : TimerMesh("m-cell-center", m_) {}
  void F() override {
    int a = 0;
    for (auto c : m.Cells()) {
      a += m.GetCenter(c)[0];
    }
    EXPOSE(a);
  }
};

class CellNCell : public TimerMesh {
 public:
  CellNCell(M& m_) : TimerMesh("m-cell-n-cell", m_) {}
  void F() override {
    size_t a = 1;
    for (auto c : m.Cells()) {
      for (auto q : m.Nci(c)) {
        a += m.GetCell(c, q).GetRaw();
      }
    }
    EXPOSE(a);
  }
};

class CellNFace : public TimerMesh {
 public:
  CellNFace(M& m_) : TimerMesh("m-cell-n-face", m_) {}
  void F() override {
    size_t a = 1;
    for (auto c : m.Cells()) {
      for (auto q : m.Nci(c)) {
        a += m.GetFace(c, q).GetRaw();
      }
    }
    EXPOSE(a);
  }
};

class CellOutward : public TimerMesh {
 public:
  CellOutward(M& m_) : TimerMesh("m-cell-outward", m_) {}
  void F() override {
    size_t a = 0;
    for (auto c : m.Cells()) {
      for (auto q : m.Nci(c)) {
        a += m.GetOutwardFactor(c, q);
      }
    }
    EXPOSE(a);
  }
};

class CellNNode : public TimerMesh {
 public:
  CellNNode(M& m_) : TimerMesh("m-cell-n-node", m_) {}
  void F() override {
    size_t a = 1;
    for (auto c : m.Cells()) {
      for (size_t q = 0; q < m.GetNumNodes(c); ++q) {
        a += m.GetNode(c, q).GetRaw();
      }
    }
    EXPOSE(a);
  }
};

class FaceCenter : public TimerMesh {
 public:
  FaceCenter(M& m_) : TimerMesh("m-face-center", m_) {}
  void F() override {
    int a = 0;
    for (auto f : m.Faces()) {
      a += m.GetCenter(f)[0];
    }
    EXPOSE(a);
  }
};

class FaceSurf : public TimerMesh {
 public:
  FaceSurf(M& m_) : TimerMesh("m-face-surf", m_) {}
  void F() override {
    int a = 0;
    for (auto f : m.Faces()) {
      a += m.GetSurface(f)[0];
    }
    EXPOSE(a);
  }
};

class FaceArea : public TimerMesh {
 public:
  FaceArea(M& m_) : TimerMesh("m-face-area", m_) {}
  void F() override {
    int a = 0;
    for (auto f : m.Faces()) {
      a += m.GetSurface(f)[0];
    }
    EXPOSE(a);
  }
};

class FaceNCell : public TimerMesh {
 public:
  FaceNCell(M& m_) : TimerMesh("m-face-n-cell", m_) {}
  void F() override {
    size_t a = 1;
    for (auto f : m.Faces()) {
      for (size_t q = 0; q < m.GetNumCells(f); ++q) {
        a += m.GetCell(f, q).GetRaw();
      }
    }
    EXPOSE(a);
  }
};

class FaceNNode : public TimerMesh {
 public:
  FaceNNode(M& m_) : TimerMesh("m-face-n-node", m_) {}
  void F() override {
    size_t a = 1;
    for (auto f : m.Faces()) {
      for (size_t q = 0; q < m.GetNumNodes(f); ++q) {
        a += m.GetNode(f, q).GetRaw();
      }
    }
    EXPOSE(a);
  }
};

class Interp : public TimerMesh {
 public:
  Interp(M& m_) : TimerMesh("interp", m_), fc(m), ff(m) {
    for (auto i : m.AllCells()) {
      fc[i] = std::sin(i.GetRaw());
    }
    auto& bf = m.GetIndexFaces();
    for (auto f : m.Faces()) {
      if (bf.GetMIdx(f)[0] == 0 && bf.GetDir(f) == Dir::i) {
        mfc[f].type = BCondType::neumann;
      }
    }
    assert(mfc.size() > 0);
  }
  void F() override {
    size_t a = 0;
    ff = UEmbed<M>::Interpolate(fc, mfc, m);
    a = ff[IdxFace(a)];
    EXPOSE(a);
  }

 private:
  FieldCell<Scal> fc;
  MapEmbed<BCond<Scal>> mfc;
  FieldFace<Scal> ff;
};

class Grad : public TimerMesh {
 public:
  Grad(M& m_) : TimerMesh("grad", m_), fc(m), ff(m) {
    for (auto i : m.AllFaces()) {
      ff[i] = std::sin(i.GetRaw());
    }
  }
  void F() override {
    static size_t a = 0;
    fc = UEmbed<M>::Gradient(ff, m);
    a += fc[IdxCell(0)][0];
    EXPOSE(a);
  }

 private:
  FieldCell<Vect> fc;
  FieldFace<Scal> ff;
};

class ExplVisc : public TimerMesh {
 public:
  ExplVisc(M& m_) : TimerMesh("explvisc", m_), fcv(m), fcf(m), ffmu(m) {
    for (auto i : m.AllCells()) {
      auto a = i.GetRaw();
      fcv[i] = Vect(std::sin(a), std::sin(a + 1), std::sin(a + 2));
    }
    for (auto i : m.AllFaces()) {
      auto a = i.GetRaw();
      ffmu[i] = std::sin(a);
    }
  }
  void F() override {
    static size_t a = 0;
    for (size_t n = 0; n < dim; ++n) {
      FieldCell<Scal> fc = GetComponent(fcv, n);
      auto ff = UEmbed<M>::Interpolate(fc, {}, m);
      auto gc = UEmbed<M>::Gradient(ff, m);
      auto gf = UEmbed<M>::Interpolate(gc, {}, m); // adhoc: zero-der cond
      for (auto idxcell : m.SuCells()) {
        Vect sum(0);
        for (size_t i = 0; i < m.GetNumFaces(idxcell); ++i) {
          IdxFace idxface = m.GetFace(idxcell, i);
          sum += gf[idxface] *
                 (ffmu[idxface] * m.GetOutwardSurface(idxcell, i)[n]);
        }
        fcf[idxcell] += sum / m.GetVolume(idxcell);
      }
    }
    a += fcf[IdxCell(0)][0];
    EXPOSE(a);
  }

 private:
  FieldCell<Vect> fcv;
  FieldCell<Vect> fcf;
  FieldFace<Scal> ffmu;
};

// i: target index
// k: current index
// Output:
// create if k == i
// ++k
// p: pointer to new instance
template <class T>
void Try(M& m, size_t i, size_t& k, Timer*& p) {
  if (k++ == i) {
    p = new T(m);
  }
}

// i: test index
// m: mesh
// Output:
// t: total per one call [sec]
// n: number of calls
// mem: memory usage in bytes
// name: test name
// Returns 1 if test with index i found
bool Run(
    const size_t i, M& m, double& t, size_t& n, size_t& mem,
    std::string& name) {
  size_t k = 0;
  Timer* p = nullptr;

  Try<Empty>(m, i, k, p);
  Try<LoopPlain>(m, i, k, p);
  Try<LoopAllCells>(m, i, k, p);
  Try<LoopInCells>(m, i, k, p);
  Try<LoopAllFaces>(m, i, k, p);
  Try<LoopInFaces>(m, i, k, p);
  Try<LoopFldPlain>(m, i, k, p);
  Try<LoopFldAllCells>(m, i, k, p);
  Try<LoopMIdxAllCells>(m, i, k, p);
  Try<LoopMIdxAllFaces>(m, i, k, p);
  Try<Interp>(m, i, k, p);
  Try<Grad>(m, i, k, p);
  Try<ExplVisc>(m, i, k, p);

  Try<CellVolume>(m, i, k, p);
  Try<CellCenter>(m, i, k, p);
  Try<FaceCenter>(m, i, k, p);
  Try<FaceSurf>(m, i, k, p);
  Try<FaceArea>(m, i, k, p);
  Try<CellNCell>(m, i, k, p);
  Try<CellNFace>(m, i, k, p);
  Try<CellOutward>(m, i, k, p);
  Try<CellNNode>(m, i, k, p);
  Try<FaceNCell>(m, i, k, p);
  Try<FaceNNode>(m, i, k, p);

  if (!p) {
    return false;
  }

  std::pair<double, size_t> e = p->Run();
  t = e.first;
  n = e.second;
  mem = sysinfo::GetMem();
  name = p->GetName();
  delete p;

  return true;
}

int main() {
  // mesh size
  std::vector<MIdx> ss = {MIdx(4),  MIdx(8),  MIdx(16),
                          MIdx(32), MIdx(64), MIdx(128)};

  using std::setw;
  std::cout << setw(20) << "name" << setw(20) << "t/allcells [ns]" << setw(20)
            << "t/incells [ns]" << setw(20) << "t [ns]" << setw(20) << "iters"
            << setw(20) << "mem [MB]" << setw(20) << "mem/allcells [B]"
            << std::endl
            << std::endl;

  for (auto s : ss) {
    size_t mem0 = sysinfo::GetMem();
    auto m = GetMesh(s);
    const size_t nca = m.GetAllBlockCells().size();
    const size_t nci = m.GetInBlockCells().size();
    std::cout << "Mesh"
              << " size=" << s << " allcells=" << nca << " incells=" << nci
              << std::endl;

    int i = 0;
    double t;
    size_t n;
    size_t mem;
    std::string name;
    while (Run(i++, m, t, n, mem, name)) {
      size_t dmem = mem - mem0;
      std::cout << setw(20) << name << setw(20) << t * 1e9 / nca << setw(20)
                << t * 1e9 / nci << setw(20) << t * 1e9 * n << setw(20) << n
                << setw(20) << (dmem / double(1 << 20)) << setw(20)
                << (dmem / nca) << std::endl;
    }
    std::cout << std::endl;
  }
}
