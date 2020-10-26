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

template <class M>
void CalcNormalHeight2(
    M& m, const FieldCell<Scal>& fcu, const FieldCell<bool>& fci, size_t edim,
    bool force_overwrite, FieldCell<Vect>& fcn) {
  fcn.Reinit(m);

  using Direction = typename M::Direction;
  for (auto c : m.SuCellsM()) {
    if (!fci[c]) {
      continue;
    }
    Vect best_n(0);
    size_t best_dz(0);
    for (size_t ddz : {0, 1, 2}) {
      if (ddz >= edim) {
        continue;
      }
      const Direction dz(ddz);
      const auto dx = dz >> 1;
      const auto dy = dz >> 2;

      auto hh = [&](Direction d) {
        return fcu[c + d - dz] + fcu[c + d] + fcu[c + d + dz];
      };

      Vect n;
      n[dx] = hh(dx) - hh(-dx);
      n[dy] = hh(dy) - hh(-dy);
      n[dz] = (fcu[c + dz] - fcu[c - dz] > 0 ? 2 : -2);
      n /= -n.norm1();
      if (std::abs(best_n[best_dz]) < std::abs(n[dz])) {
        best_n = n;
        best_dz = dz;
      }
    }

    if (force_overwrite ||
        std::abs(best_n[best_dz]) < std::abs(fcn[c][best_dz])) {
      fcn[c] = best_n;
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

#define CMP(fca, fcb)                                                      \
  do {                                                                     \
    Scal diff;                                                             \
    if ((diff = DiffMax(fca, fcb, m)) > eps) {                             \
      std::cerr << util::Format(                                           \
          "Verify {}: max difference ({},{}) exceeded: {:.4e} > {:.4e}\n", \
          name, #fca, #fcb, diff, eps);                                    \
      ++failed;                                                            \
    }                                                                      \
  } while (0);

// Returns the number of failed tests
template <int test>
static int Verify() {
  int failed = 0;
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
    CalcNormalHeight2(m, fcu, fcmask, edim, true, fcn2);
    CMP(fcn, fcn1);
    CMP(fcn, fcn2);
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
    sizes.emplace_back(n, n, n);
  }

  const std::string fmt = "{:10} {:16} {:16}\n";

  for (auto size : sizes) {
    auto m = GetMesh(size);
    const size_t sucells = m.GetSuBlockCells().size();
    std::cout << "Mesh"
              << " size=" << size << " sucells=" << sucells << std::endl;

    int test = 0;
    double time;
    size_t niters;
    size_t mem;
    std::string name;
    std::cout << util::Format(fmt, "name", "time/sucells", "niters");
    while (Run(test++, m, time, niters, mem, name)) {
      std::cout << util::Format(fmt, name, time * 1e9 / sucells, niters);
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }
  if (failed) {
    std::cerr << util::Format("{} tests failed", failed) << std::endl;
  }
  return failed;
}
