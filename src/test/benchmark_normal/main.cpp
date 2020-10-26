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
typename M::Scal DiffMax(
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

Scal Random(Scal q) {
  return std::sin(std::sin(q * 123.456) * 654.321);
}

Scal RandomNext(Scal& q) {
  return Random(q += 0.1);
}

Mesh GetMesh(MIdx size) {
  const Rect<Vect> dom(Vect(0), Vect(1));
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
    } else {
      Normal::CalcNormalHeight1(m, fcu, fcmask, edim, true, fcn);
    }

    ii = fcn[IdxCell(ii)][0];
  }

 private:
  FieldCell<Scal> fcu;
  FieldCell<bool> fcmask;
  FieldCell<Vect> fcn;
  FieldCell<char> fcd;
};

class Partstr : public TimerMesh {
 public:
  Partstr(Mesh& m_)
      : TimerMesh("partstr", m_), fcu(m), fcmask(m, true), fcn(m), fca(m) {
    for (auto c : m.AllCells()) {
      Scal q = c.GetRaw();
      fcu[c] = RandomNext(q);
      fcn[c] = Vect(RandomNext(q), RandomNext(q), RandomNext(q));
      fca[c] = RandomNext(q);
    }
  }
  void F() override {
    volatile size_t ii = 0;

    ii = fcn[IdxCell(ii)][0];
  }

 private:
  FieldCell<Scal> fcu;
  FieldCell<bool> fcmask;
  FieldCell<Vect> fcn;
  FieldCell<Scal> fca;
};

// f=0: youngs
// f=1: height edim=3
// f=2: height edim=2
void Cmp(int f) {
  auto m = GetMesh(MIdx(8));
  FieldCell<Scal> fcu(m);
  FieldCell<bool> fcmask(m, true);
  FieldCell<Vect> fcn(m);
  FieldCell<Vect> fcn2(m);
  for (auto c : m.AllCells()) {
    fcu[c] = Random(c.GetRaw());
  }
  if (f == 0) {
    Normal::CalcNormalYoung(m, fcu, fcmask, fcn);
    Normal::CalcNormalYoung1(m, fcu, fcmask, fcn2);
  } else if (f == 1) {
    size_t edim = 3;
    Normal::CalcNormalHeight(m, fcu, fcmask, edim, true, fcn);
    Normal::CalcNormalHeight1(m, fcu, fcmask, edim, true, fcn2);
  } else {
    size_t edim = 2;
    Normal::CalcNormalHeight(m, fcu, fcmask, edim, true, fcn);
    Normal::CalcNormalHeight1(m, fcu, fcmask, edim, true, fcn2);
  }

  const Scal r = DiffMax(fcn, fcn2, m);
  std::vector<std::string> names = {
      "NormalYoung", "NormalHeight,edim=3", "NormalHeight,edim=2"};
  const Scal eps = 1e-15;
  fassert(
      r <= eps, //
      util::Format(
          "Cmp {}: max difference exceeded: {:.16e} > {:.16e}", names[f], r,
          eps));
}

// test: test index
// Output:
// time: total per one call [sec]
// niters: number of calls
// mem: memory usage in bytes
// name: test name
// Returns 1 if test with index i found
bool Run(
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
  create((Partstr*)0);

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
  Cmp(0);
  Cmp(1);
  Cmp(2);

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
}
