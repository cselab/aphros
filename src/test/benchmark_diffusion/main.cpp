// Created by Petr Karnakov on 26.04.2019
// Copyright 2019 ETH Zurich

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

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

M GetMesh(MIdx size) {
  Rect<Vect> dom(Vect(0), Vect(1));
  MIdx begin(0, 0, 0); // lower index
  int halos = 2;
  return {begin, size, dom, halos, true, true, size, 0};
}

class TimerMesh : public ExecutionTimer {
 public:
  TimerMesh(const std::string& name, M& m_)
      : ExecutionTimer(name, 0.1, 3), m(m_) {}

 protected:
  M& m;
};

class LoopPlain : public TimerMesh {
 public:
  LoopPlain(M& m_) : TimerMesh("loop-plain", m_) {}
  void F() override {
    volatile size_t a = 0;
    size_t b = a;
    for (size_t i = 0, ie = m.GetInBlockCells().size(); i < ie; ++i) {
      b += i;
    }
    a = b;
  }
};

class LoopInCells : public TimerMesh {
 public:
  LoopInCells(M& m_) : TimerMesh("loop-incells", m_) {}
  void F() override {
    volatile size_t a = 0;
    size_t b = a;
    for (auto i : m.Cells()) {
      b += i.GetRaw();
    }
    a = b;
  }
};

class LoopFldInCells : public TimerMesh {
 public:
  LoopFldInCells(M& m_) : TimerMesh("loop-fld-incells", m_), v(m) {}
  void F() override {
    volatile size_t a = 0;
    size_t b = a;
    for (auto i : m.AllCells()) {
      b += v[i];
    }
    a = b;
    v[IdxCell(a % v.size())] = a;
  }

 private:
  FieldCell<Scal> v;
};

class Interp : public TimerMesh {
 public:
  Interp(M& m_) : TimerMesh("interp", m_), fc(m), ff(m) {
    using Dir = typename M::Dir;
    for (auto i : m.AllCells()) {
      fc[i] = std::sin(i.GetRaw());
    }
    auto& bf = m.GetIndexFaces();
    for (auto f : m.Faces()) {
      if (bf.GetMIdx(f)[0] == 0 && bf.GetDir(f) == Dir(0)) {
        mfc[f].type = BCondType::neumann;
      }
    }
    assert(mfc.size() > 0);
  }
  void F() override {
    volatile size_t a = 0;
    ff = UEmbed<M>::Interpolate(fc, mfc, m);
    a = ff[IdxFace(a)];
  }

 private:
  FieldCell<Scal> fc;
  MapEmbed<BCond<Scal>> mfc;
  FieldFace<Scal> ff;
};

class Diffusion : public TimerMesh {
 public:
  Diffusion(M& m_) : TimerMesh("diff-face-idx", m_), fc(m), ff(m) {
    for (auto i : m.AllCells()) {
      fc[i] = std::sin(i.GetRaw());
    }
  }
  void F() override {
    volatile size_t ii = 0;

    // normal gradient
    for (auto f : m.Faces()) {
      auto cm = m.GetCell(f, 0);
      auto cp = m.GetCell(f, 1);
      ff[f] = (fc[cp] - fc[cm]);
    }

    // advance
    for (auto c : m.Cells()) {
      Scal s = 0;
      for (auto q : m.Nci(c)) {
        IdxFace f = m.GetFace(c, q);
        s += ff[f] * m.GetOutwardFactor(c, q);
      }
      fc[c] += s;
    }

    ii = fc[IdxCell(ii)];
  }

 private:
  FieldCell<Scal> fc;
  FieldFace<Scal> ff;
};

class DiffusionNeighb : public TimerMesh {
 public:
  DiffusionNeighb(M& m_) : TimerMesh("diff-cell-idx", m_), fc(m), fct(m) {
    for (auto i : m.AllCells()) {
      fc[i] = std::sin(i.GetRaw());
    }
    fct = fc;
  }
  void F() override {
    volatile size_t ii = 0;

    fc.swap(fct);
    for (auto c : m.Cells()) {
      Scal s = 0;
      for (auto q : m.Nci(c)) {
        IdxCell cc = m.GetCell(c, q);
        s += fct[cc];
      }
      s += -6 * fct[c];
      fc[c] = fct[c] + s;
    }

    ii = fc[IdxCell(ii)];
  }

 private:
  FieldCell<Scal> fc;
  FieldCell<Scal> fct;
};

class DiffusionPlain : public TimerMesh {
 public:
  DiffusionPlain(M& m_) : TimerMesh("diff-cell-plain", m_), fc(m), fct(m) {
    for (auto i : m.AllCells()) {
      fc[i] = std::sin(i.GetRaw());
    }
    fct = fc;
  }
  void F() override {
    volatile size_t ii = 0;

    auto& bc = m.GetIndexCells();
    MIdx wb = bc.GetBegin();
    MIdx ws = bc.GetSize();
    const size_t hl = -wb[0]; // halo cells
    const size_t nx = ws[0]; // block size n^3
    const size_t ny = ws[1]; // block size n^3
    const size_t nz = ws[2]; // block size n^3

    fc.swap(fct);
    Scal* d = fc.data();
    Scal* t = fct.data();
    const size_t dx = 1;
    const size_t dy = nx;
    const size_t dz = nx * ny;
    for (size_t z = hl, ze = nz - hl; z < ze; ++z) {
      for (size_t y = hl, ye = ny - hl; y < ye; ++y) {
        size_t i = (z * ny + y) * nx + hl;
        for (size_t x = hl, xe = nx - hl; x < xe; ++x) {
          Scal s = -6 * t[i];
          s += t[i + dx];
          s += t[i - dx];
          s += t[i + dy];
          s += t[i - dy];
          s += t[i + dz];
          s += t[i - dz];
          d[i] = t[i] + s;
          ++i;
        }
      }
    }

    ii = fc[IdxCell(ii)];
  }

 private:
  FieldCell<Scal> fc;
  FieldCell<Scal> fct;
};

class DiffusionPlainFace : public TimerMesh {
 public:
  DiffusionPlainFace(M& m_)
      : TimerMesh("diff-face-plain", m_), fc(m), fcx(m), fcy(m), fcz(m) {
    for (auto i : m.AllCells()) {
      fc[i] = std::sin(i.GetRaw());
    }
  }
  void F() override {
    volatile size_t ii = 0;

    auto& bc = m.GetIndexCells();
    MIdx wb = bc.GetBegin();
    MIdx ws = bc.GetSize();
    const size_t hl = -wb[0]; // halo cells
    const size_t nx = ws[0]; // block size n^3
    const size_t ny = ws[1]; // block size n^3
    const size_t nz = ws[2]; // block size n^3

    Scal* d = fc.data();
    Scal* tx = fcx.data();
    Scal* ty = fcy.data();
    Scal* tz = fcz.data();

    const size_t dx = 1;
    const size_t dy = nx;
    const size_t dz = nx * ny;

    for (size_t z = hl, ze = nz - hl; z < ze; ++z) {
      for (size_t y = hl, ye = ny - hl; y < ye; ++y) {
        for (size_t x = hl, xe = nx - hl; x <= xe; ++x) {
          size_t i = (z * ny + y) * nx + x;
          tx[i] = d[i + dx] - d[i];
        }
      }
    }
    for (size_t z = hl, ze = nz - hl; z < ze; ++z) {
      for (size_t y = hl, ye = ny - hl; y <= ye; ++y) {
        for (size_t x = hl, xe = nx - hl; x < xe; ++x) {
          size_t i = (z * ny + y) * nx + x;
          ty[i] = d[i + dy] - d[i];
        }
      }
    }
    for (size_t z = hl, ze = nz - hl; z <= ze; ++z) {
      for (size_t y = hl, ye = ny - hl; y < ye; ++y) {
        for (size_t x = hl, xe = nx - hl; x < xe; ++x) {
          size_t i = (z * ny + y) * nx + x;
          tz[i] = d[i + dz] - d[i];
        }
      }
    }
    for (size_t z = hl, ze = nz - hl; z < ze; ++z) {
      for (size_t y = hl, ye = ny - hl; y < ye; ++y) {
        for (size_t x = hl, xe = nx - hl; x < xe; ++x) {
          size_t i = (z * ny + y) * nx + x;
          d[i] = tx[i + dx] - tx[i] + ty[i + dy] - ty[i] + tz[i + dz] - tz[i];
        }
      }
    }

    ii = fc[IdxCell(ii)];
  }

 private:
  FieldCell<Scal> fc;
  FieldCell<Scal> fcx;
  FieldCell<Scal> fcy;
  FieldCell<Scal> fcz;
};

// i: target index
// k: current index
// Output:
// create if k == i
// ++k
// p: pointer to new instance
template <class T>
void Try(M& m, size_t i, size_t& k, ExecutionTimer*& p) {
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
  ExecutionTimer* p = nullptr;

  // Try<LoopPlain>(m, i, k, p);
  // Try<LoopInCells>(m, i, k, p);
  // Try<LoopFldInCells>(m, i, k, p);
  // Try<Interp>(m, i, k, p);
  Try<Diffusion>(m, i, k, p);
  Try<DiffusionPlainFace>(m, i, k, p);
  Try<DiffusionNeighb>(m, i, k, p);
  Try<DiffusionPlain>(m, i, k, p);

  if (!p) {
    return false;
  }

  auto e = p->Run();
  t = e.min_call_time;
  n = e.iters;
  mem = sysinfo::GetMem();
  name = p->GetName();
  delete p;

  return true;
}

int main() {
  // mesh size
  std::vector<MIdx> ss;
  for (int n : {16, 32, 64, 128}) {
    ss.emplace_back(n, n, n);
  }

  const size_t ww = 16;

  using std::setw;

  for (auto s : ss) {
    auto m = GetMesh(s);
    const size_t nci = m.GetInBlockCells().size();
    std::cout << "Mesh"
              << " size=" << s << " incells=" << nci << std::endl;

    /*
    std::cout
        << setw(ww) << "name"
        << setw(ww) << "t/incells [ns]"
        << setw(ww) << "iters"
        << std::endl;
    for (size_t q = 0; q < ww * 3; ++q) {
      std::cout << "-";
    }
    std::cout << std::endl;
    */

    int i = 0;
    double t;
    size_t n;
    size_t mem;
    std::string name;
    while (Run(i++, m, t, n, mem, name)) {
      std::cout << setw(ww) << name << setw(ww) << t * 1e9 / nci << setw(ww)
                << n << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }
}
