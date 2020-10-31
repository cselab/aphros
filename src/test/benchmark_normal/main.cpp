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
    if (IsNan(u[i] - v[i])) {
      return GetNan<Scal>();
    }
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

#include <x86intrin.h>

std::ostream& operator<<(std::ostream& out, const __m256d& d) {
  constexpr int width = 4;
  double dd[width];
  _mm256_storeu_pd(dd, d);
  bool first = true;
  for (size_t i = 0; i < width; ++i) {
    if (!first) {
      out << " ";
    }
    first = false;
    out << dd[i];
  }
  return out;
}


template <class M>
static void CalcNormalYoung2(
    M& m, const FieldCell<Scal>& fcu, const FieldCell<bool>& fci,
    FieldCell<Vect>& fcn) {
  using MIdx = typename M::MIdx;
  auto ic = m.GetIndexCells();
  auto bc = m.GetSuBlockCells();
  const MIdx s = ic.GetSize();
  const size_t nx = s[0];
  const size_t ny = s[1];
  // offset
  const size_t fx = 1;
  const size_t fy = nx;
  const size_t fz = ny * nx;

  fcn.Reinit(m, Vect(0));

  // index range
  const MIdx wb = bc.GetBegin() - ic.GetBegin();
  const MIdx we = bc.GetEnd() - ic.GetBegin();
  const Scal* pu = fcu.data();
  Vect* pn = fcn.data();
  const bool* pi = fci.data();

  auto interleave = [](const __m256d& x, const __m256d& y, const __m256d& z,
                       __m256d& v0, __m256d& v1, __m256d& v2) {
    // x: x0 x1 x2 x3
    // y: y0 y1 y2 y3
    // z: z0 z1 z2 z3
    auto xp = _mm256_permute2f128_pd(x, x, 0b00000001); // x2 x3 x0 x1
    auto zp = _mm256_permute2f128_pd(z, z, 0b00000001); // z2 z3 z0 z1
    auto xs = _mm256_blend_pd(x, xp, 0b1010); /// x0 x3 x2 x1
    auto ys = _mm256_shuffle_pd(y, y, 0b0101); // y1 y0 y3 y2
    auto zs = _mm256_blend_pd(z, zp, 0b0101); /// z2 z1 z0 z3
    auto xyb = _mm256_blend_pd(xs, ys, 0b1010); // x0 y0 x2 y2
    auto yzb = _mm256_blend_pd(ys, zs, 0b1010); // y1 z1 y3 z3
    auto zxb = _mm256_blend_pd(zs, xs, 0b1010); // z2 x3 z0 x1
    v0 = _mm256_blend_pd(xyb, zxb, 0b1100); // x0 y0 z0 x1
    v1 = _mm256_blend_pd(yzb, xyb, 0b1100); // y1 z1 x2 y2
    v2 = _mm256_blend_pd(zxb, yzb, 0b1100); // y2 z3 x3 y3
  };

  auto normalize1 = [](__m256d& x, __m256d& y, __m256d& z) {
    auto xa = _mm256_max_pd(_mm256_sub_pd(_mm256_setzero_pd(), x), x);
    auto ya = _mm256_max_pd(_mm256_sub_pd(_mm256_setzero_pd(), y), y);
    auto za = _mm256_max_pd(_mm256_sub_pd(_mm256_setzero_pd(), z), z);
    auto sum = _mm256_add_pd(xa, ya);
    sum = _mm256_add_pd(sum, za);
    sum = _mm256_sub_pd(_mm256_setzero_pd(), sum);
    x = _mm256_div_pd(x, sum);
    y = _mm256_div_pd(y, sum);
    z = _mm256_div_pd(z, sum);
  };

  const __m256d c1 = _mm256_set1_pd(1);
  const __m256d c2 = _mm256_set1_pd(2);
  const __m256d c4 = _mm256_set1_pd(4);
  constexpr int W = 4;
  for (int z = wb[2]; z < we[2]; ++z) {
    for (int y = wb[1]; y < we[1]; ++y) {
      for (int x = wb[0]; x < we[0]; x += W) {
        int i = (z * ny + y) * nx + x;
        if (!pi[i]) {
          continue;
        }
        auto q = [i, fy, fz, pu](int dx, int dy, int dz) -> const Scal* {
          size_t ii = i + dx * fx + dy * fy + dz * fz;
          return pu + ii;
        };
        auto add = [&q](
                       __m256d& vx, __m256d& vy, __m256d& vz, const __m256d k,
                       int dx, int dy, int dz) {
          const __m256d mqx = _mm256_loadu_pd(q(dx, dy, dz));
          const __m256d mqy = _mm256_loadu_pd(q(dz, dx, dy));
          const __m256d mqz = _mm256_loadu_pd(q(dy, dz, dx));
          vx = _mm256_fmadd_pd(k, mqx, vx);
          vy = _mm256_fmadd_pd(k, mqy, vy);
          vz = _mm256_fmadd_pd(k, mqz, vz);
        };
        auto sub = [&q](
                       __m256d& vx, __m256d& vy, __m256d& vz, const __m256d k,
                       int dx, int dy, int dz) {
          const __m256d mqx = _mm256_loadu_pd(q(dx, dy, dz));
          const __m256d mqy = _mm256_loadu_pd(q(dz, dx, dy));
          const __m256d mqz = _mm256_loadu_pd(q(dy, dz, dx));
          vx = _mm256_fnmadd_pd(k, mqx, vx);
          vy = _mm256_fnmadd_pd(k, mqy, vy);
          vz = _mm256_fnmadd_pd(k, mqz, vz);
        };
        __m256d vx = _mm256_setzero_pd();
        __m256d vy = _mm256_setzero_pd();
        __m256d vz = _mm256_setzero_pd();
        sub(vx, vy, vz, c1, -1, -1, -1);
        sub(vx, vy, vz, c2, -1, -1, +0);
        sub(vx, vy, vz, c1, -1, -1, +1);
        sub(vx, vy, vz, c2, -1, +0, -1);
        sub(vx, vy, vz, c4, -1, +0, +0);
        sub(vx, vy, vz, c2, -1, +0, +1);
        sub(vx, vy, vz, c1, -1, +1, -1);
        sub(vx, vy, vz, c2, -1, +1, +0);
        sub(vx, vy, vz, c1, -1, +1, +1);
        add(vx, vy, vz, c1, +1, -1, -1);
        add(vx, vy, vz, c2, +1, -1, +0);
        add(vx, vy, vz, c1, +1, -1, +1);
        add(vx, vy, vz, c2, +1, +0, -1);
        add(vx, vy, vz, c4, +1, +0, +0);
        add(vx, vy, vz, c2, +1, +0, +1);
        add(vx, vy, vz, c1, +1, +1, -1);
        add(vx, vy, vz, c2, +1, +1, +0);
        add(vx, vy, vz, c1, +1, +1, +1);

        normalize1(vx, vy, vz);

        __m256d d0;
        __m256d d1;
        __m256d d2;

        interleave(vx, vy, vz, d0, d1, d2);
        _mm256_storeu_pd((Scal*)(pn + i) + 0, d0);
        _mm256_storeu_pd((Scal*)(pn + i) + 4, d1);
        _mm256_storeu_pd((Scal*)(pn + i) + 8, d2);
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
class Youngs : public TimerMesh {
 public:
  Youngs(Mesh& m_)
      : TimerMesh("young" + std::to_string(id), m_), fcu(m), fcmask(m, true) {
    for (auto i : m.AllCells()) {
      fcu[i] = Random(i.GetRaw());
    }
  }
  void F() override {
    volatile size_t ii = 0;

    if (id == 0) {
      Normal::CalcNormalYoung(m, fcu, fcmask, fcn);
    } else if (id == 1) {
      Normal::CalcNormalYoung1(m, fcu, fcmask, fcn);
    } else if (id == 2) {
      CalcNormalYoung2(m, fcu, fcmask, fcn);
    } else {
      fassert(false);
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
    } else if (id == 2) {
      CalcNormalHeight2(m, fcu, fcmask, edim, true, fcn);
    } else {
      fassert(false);
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
    const Scal diff = DiffMax(fca, fcb, m);                                \
    if (!(diff <= eps)) {                                                  \
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
    CalcNormalYoung2(m, fcu, fcmask, fcn2);
    CMP(fcn, fcn1);
    CMP(fcn, fcn2);
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

  create((Youngs<0>*)0);
  create((Youngs<1>*)0);
  create((Youngs<2>*)0);
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

  const std::string fmt = "{:10} {:16.2f} {:16}\n";

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
    std::cout << util::Format(fmt, "name", "ns/cell", "niters");
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
