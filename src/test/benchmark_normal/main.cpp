// Created by Petr Karnakov on 30.04.2019
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
#include "solver/normal.h"
#include "solver/normal.ipp"
#include "solver/solver.h"
#include "util/format.h"
#include "util/sysinfo.h"
#include "util/timer.h"

template <class Idx, class M>
static typename M::Scal DiffMax(
    const GField<typename M::Vect, Idx>& u,
    const GField<typename M::Vect, Idx>& v, const M& m) {
  using Scal = typename M::Scal;
  Scal r = 0;
  for (auto i : m.template GetIn<Idx>()) {
    r = std::max(r, (u[i] - v[i]).norminf());
  }
  return r;
}

const int dim = 3;
using MIdx = GMIdx<dim>;
using IdxCell = IdxCell;
using IdxFace = IdxFace;
using Dir = GDir<dim>;
using Scal = double;
using Vect = generic::Vect<Scal, dim>;
using Mesh = MeshStructured<Scal, dim>;
using Normal = typename UNormal<Mesh>::Imp;

template <class Q>
static void F(
    Q& qn, Scal& gxm, Scal& gxp, Scal& gym, Scal& gyp, Scal& sg, Scal& gc,
    Scal& gmm, Scal& gmp, Scal& gpm, Scal& gpp) {
  auto g = [&qn](int dx, int dy) {
    return qn(dx, dy, 1) + qn(dx, dy, 0) + qn(dx, dy, -1);
  };
  gxm = g(-1, 0);
  gxp = g(1, 0);
  gym = g(0, -1);
  gyp = g(0, 1);
  sg = (qn(0, 0, 1) - qn(0, 0, -1) > 0. ? 1. : -1);
  gc = g(0, 0);
  gmm = g(-1, -1);
  gmp = g(-1, 1);
  gpm = g(1, -1);
  gpp = g(1, 1);
}

template <class Q>
// TODO: check what optimizer generates with lambdas
static void GetHeight(Q& q, Vect& pn, size_t edim, bool ow) {
  auto qx = [&q](int dx, int dy, int dz) { return q(dz, dx, dy); };
  auto qy = [&q](int dx, int dy, int dz) { return q(dy, dz, dx); };
  auto qz = [&q](int dx, int dy, int dz) { return q(dx, dy, dz); };

  Vect bn(0); // best normal
  size_t bd = 0; // best direction
  std::vector<size_t> dd; // direction of plane normal
  if (edim == 2) {
    dd = {0, 1};
  } else {
    dd = {0, 1, 2};
  }
  for (size_t dn : dd) {
    size_t dx = (dn + 1) % dim;
    size_t dy = (dn + 2) % dim;

    // height function
    Scal gxm, gxp, gym, gyp;
    Scal gc, gmm, gmp, gpm, gpp;
    // sign: +1 if u increases in dn
    Scal sg;
    switch (dn) {
      case 0: { // dx:1 , dy:2
        F(qx, gxm, gxp, gym, gyp, sg, gc, gmm, gmp, gpm, gpp);
        break;
      }
      case 1: { // dx:2 , dy:0
        F(qy, gxm, gxp, gym, gyp, sg, gc, gmm, gmp, gpm, gpp);
        break;
      }
      default: { // dx:0 , dy:1
        F(qz, gxm, gxp, gym, gyp, sg, gc, gmm, gmp, gpm, gpp);
        break;
      }
    }

    // first derivative (slope)
    Scal gx = (gxp - gxm) * 0.5; // centered
    Scal gy = (gyp - gym) * 0.5;
    // outer normal
    Vect n;
    n[dx] = -gx;
    n[dy] = -gy;
    n[dn] = -sg;
    n /= n.norm1();
    // select best with minimal slope
    if (std::abs(n[dn]) > std::abs(bn[bd])) {
      bn = n;
      bd = dn;
    }
  }

  // update if ow=1 or gives steeper profile in plane dn
  if (ow || std::abs(bn[bd]) < std::abs(pn[bd])) {
    pn = bn;
  }
}

template <class M>
void CalcNormalHeight2(
    M& m, const FieldCell<Scal>& fcu, const FieldCell<bool>& fci, size_t edim,
    bool ow, FieldCell<Vect>& fcn) {
  using MIdx = typename M::MIdx;
  auto ic = m.GetIndexCells();
  auto bc = m.GetSuBlockCells();
  MIdx s = ic.GetSize();
  const size_t nx = s[0];
  const size_t ny = s[1];
  // offset
  const size_t fx = 1;
  const size_t fy = nx;
  const size_t fz = ny * nx;

  // index range
  const MIdx wb = bc.GetBegin() - ic.GetBegin();
  const MIdx we = bc.GetEnd() - ic.GetBegin();
  const size_t xb = wb[0], yb = wb[1], zb = wb[2];
  const size_t xe = we[0], ye = we[1], ze = we[2];

  fcn.Reinit(m);

  const Scal* pu = fcu.data();
  Vect* pn = fcn.data();
  const bool* pi = fci.data();
  for (size_t z = zb; z < ze; ++z) {
    for (size_t y = yb; y < ye; ++y) {
      for (size_t x = xb; x < xe; ++x) {
        size_t i = (z * ny + y) * nx + x;
        if (!pi[i]) {
          continue;
        }
        auto q = [i, fy, fz, pu](int dx, int dy, int dz) {
          size_t ii = i + dx * fx + dy * fy + dz * fz;
          return pu[ii];
        };
        GetHeight(q, pn[i], edim, ow);
      }
    }
  }
}


Scal Random(Scal q) {
  return std::sin(std::sin(q * 123.456) * 654.321);
}

Scal RandomNext(Scal& q) {
  return Random(q += 0.1);
}

Mesh GetMesh(MIdx size) {
  const Rect<Vect> dom(Vect(0), Vect(size) / Vect(size.max()));
  const MIdx begin(0, 0, 0);
  const int halos = 2;
  return InitUniformMesh<Mesh>(dom, begin, size, halos, true, true, size, 0);
}

class TimerMesh : public ExecutionTimer {
 public:
  TimerMesh(const std::string& name, Mesh& m_)
      : ExecutionTimer(name, 0.1, 3), m(m_) {}

 protected:
  Mesh& m;
};

template <int id>
class Young : public TimerMesh {
 public:
  Young(Mesh& m_)
      : TimerMesh("young" + std::to_string(id), m_), fcu(m), fcmask(m, true) {
    for (auto i : m.AllCells()) {
      fcu[i] = Random(i.GetRaw());
    }
  }
  void F() override {
    volatile size_t ii = 0;

    if (id == 0) {
      Normal::CalcNormalYoung(m, fcu, fcmask, fcn);
    } else {
      Normal::CalcNormalYoung1(m, fcu, fcmask, fcn);
    }

    ii = fcn[IdxCell(ii)][0];
  }

 private:
  FieldCell<Scal> fcu;
  FieldCell<bool> fcmask;
  FieldCell<Vect> fcn;
};

template <int id>
class Height : public TimerMesh {
 public:
  Height(Mesh& m_)
      : TimerMesh("height" + std::to_string(id), m_), fcu(m), fcmask(m, true) {
    for (auto i : m.AllCells()) {
      fcu[i] = Random(i.GetRaw());
    }
  }
  void F() override {
    volatile size_t ii = 0;
    size_t edim = 3;

    if (id == 0) {
      Normal::CalcNormalHeight(m, fcu, fcmask, edim, true, fcn);
    } else if (id == 1) {
      Normal::CalcNormalHeight1(m, fcu, fcmask, edim, true, fcn);
    } else {
      CalcNormalHeight2(m, fcu, fcmask, edim, true, fcn);
    }

    ii = fcn[IdxCell(ii)][0];
  }

 private:
  FieldCell<Scal> fcu;
  FieldCell<bool> fcmask;
  FieldCell<Vect> fcn;
  FieldCell<char> fcd;
};

#define CMP(fca, fcb)                                                        \
  do {                                                                       \
    Scal diff;                                                               \
    if ((diff = DiffMax(fca, fcb, m)) > eps) {                               \
      std::cerr << util::Format(                                             \
          "Verify {}: max difference exceeded: {:.8e} > {:.8e}", name, diff, \
          eps);                                                              \
      ++failed;                                                              \
    }                                                                        \
  } while (0);

// Returns the number of failed tests
template <int test>
static int Verify() {
  int failed =0;
  auto m = GetMesh(MIdx(32, 16, 8));
  FieldCell<Scal> fcu(m);
  FieldCell<bool> fcmask(m, true);
  FieldCell<Vect> fcn(m);
  FieldCell<Vect> fcn1(m);
  FieldCell<Vect> fcn2(m);
  for (auto c : m.AllCells()) {
    fcu[c] = Random(c.GetRaw());
  }

  std::string name;
  const Scal eps = 1e-14;

  if (test == 0) {
    name = "NormalYoung";
    Normal::CalcNormalYoung(m, fcu, fcmask, fcn);
    Normal::CalcNormalYoung1(m, fcu, fcmask, fcn1);
    CMP(fcn, fcn1);
  } else if (test == 1) {
    name = "NormalHeight,edim=3";
    size_t edim = 3;
    Normal::CalcNormalHeight(m, fcu, fcmask, edim, true, fcn);
    Normal::CalcNormalHeight1(m, fcu, fcmask, edim, true, fcn1);
    CalcNormalHeight2(m, fcu, fcmask, edim, true, fcn2);
    CMP(fcn, fcn1);
    CMP(fcn, fcn2);
  } else if (test == 2) {
    name = "NormalHeight,edim=2";
    size_t edim = 2;
    Normal::CalcNormalHeight(m, fcu, fcmask, edim, true, fcn);
    Normal::CalcNormalHeight1(m, fcu, fcmask, edim, true, fcn1);
    CMP(fcn, fcn1);
  } else {
    fassert(false, util::Format("Unknown test={}", test));
  }

  return failed;
}

// test: test index
// Output:
// time: total per one call [sec]
// niters: number of calls
// mem: memory usage in bytes
// name: test name
// Returns 1 if test with index i found
static bool Run(
    const size_t test, Mesh& m, /*out*/ double& time, size_t& niters,
    size_t& mem, std::string& name) {
  size_t current = 0;
  std::unique_ptr<ExecutionTimer> ptr;

  auto create = [&](auto* kernel) {
    if (current++ == test) {
      using T = typename std::remove_pointer<decltype(kernel)>::type;
      ptr = std::make_unique<T>(m);
    }
  };

  create((Young<0>*)0);
  create((Young<1>*)0);
  create((Height<0>*)0);
  create((Height<1>*)0);
  create((Height<2>*)0);

  if (!ptr) {
    return false;
  }

  auto e = ptr->Run();
  time = e.min_call_time;
  niters = e.iters;
  mem = sysinfo::GetMem();
  name = ptr->GetName();

  return true;
}

int main() {
  int failed = 0;
  failed += Verify<0>();
  failed += Verify<1>();
  failed += Verify<2>();

  std::vector<MIdx> sizes;
  for (int n : {8, 16, 32, 64}) {
    sizes.emplace_back(n, n, 8);
  }

  const std::string fmt = "{:10} {:16} {:16}\n";

  for (auto size : sizes) {
    auto m = GetMesh(size);
    const size_t ncells = m.GetInBlockCells().size();
    std::cout << "Mesh"
              << " size=" << size << " incells=" << ncells << std::endl;

    int test = 0;
    double time;
    size_t niters;
    size_t mem;
    std::string name;
    std::cout << util::Format(fmt, "name", "time/ncells", "niters");
    while (Run(test++, m, time, niters, mem, name)) {
      std::cout << util::Format(fmt, name, time * 1e9 / ncells, niters);
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }
  if (failed) {
    std::cerr << util::Format("{} tests failed", failed) << std::endl;
  }
  return failed;
}
