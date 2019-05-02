#include <sstream>
#include <iostream>
#include <cassert>
#include <functional>
#include <cmath>
#include <string>
#include <memory>
#include <iomanip>

#include "geom/mesh.h"
#include "solver/solver.h"
#include "util/timer.h"
#include "util/sysinfo.h"
#include "solver/normal.h"
#include "solver/normal.ipp"

const int dim = 3;
using MIdx = GMIdx<dim>;
using IdxCell = IdxCell;
using IdxFace = IdxFace;
using Dir = GDir<dim>;
using Scal = double;
using Vect = GVect<Scal, dim>;
using Mesh = MeshStructured<Scal, dim>;
using Normal = typename solver::UNormal<Mesh>::Imp;

Scal Rnd(Scal q) {
  return std::sin(std::sin(q * 123.456) * 654.321);
}

// rnd next
Scal Rndn(Scal& q) {
  q += 0.1;
  return std::sin(std::sin(q * 123.456) * 654.321);
}


template <class Idx, class M>
typename M::Scal DiffMax(
    const GField<typename M::Vect, Idx>& u,
    const GField<typename M::Vect, Idx>& v,
    const M& m) {
  using Scal = typename M::Scal;
  Scal r = 0;
  for (auto i : m.template GetIn<Idx>()) {
     r = std::max(r, (u[i] - v[i]).norminf());
  }
  return r;
}

Mesh GetMesh(MIdx s /*size in cells*/) {
  Rect<Vect> dom(Vect(0), Vect(1));
  MIdx b(0, 0, 0); // lower index
  int hl = 2;         // halos 
  return InitUniformMesh<Mesh>(dom, b, s, hl, true, true, s, 0);
}

class TimerMesh : public Timer {
 public:
  TimerMesh(const std::string& name, Mesh& m) : Timer(name, 0.1, 3), m(m) {}

 protected:
  Mesh& m;
};

template <int id>
class Young : public TimerMesh {
 public:
  Young(Mesh& m) 
      : TimerMesh("young" + std::to_string(id), m)
      , fc(m), fci(m, true)
  {
    for (auto i : m.AllCells()) {
      fc[i] = Rnd(i.GetRaw());
    }
  }
  void F() override {
    volatile size_t ii = 0;

    if (id == 0) {
      Normal::CalcNormalYoung(m, fc, fci, fcn);
    } else {
      Normal::CalcNormalYoung1(m, fc, fci, fcn);
    }

    ii = fcn[IdxCell(ii)][0];
  }

 private:
  FieldCell<Scal> fc;
  FieldCell<bool> fci;
  FieldCell<Vect> fcn;
};

template <int id>
class Height : public TimerMesh {
 public:
  Height(Mesh& m) 
      : TimerMesh("height" + std::to_string(id), m)
      , fc(m), fci(m, true)
  {
    for (auto i : m.AllCells()) {
      fc[i] = Rnd(i.GetRaw());
    }
  }
  void F() override {
    volatile size_t ii = 0;
    size_t edim = 3;

    if (id == 0) {
      Normal::CalcNormalHeight(m, fc, fci, edim, true, fcn, fck);
    } else {
      Normal::CalcNormalHeight1(m, fc, fci, edim, true, fcn, fck);
    }

    ii = fcn[IdxCell(ii)][0];
  }

 private:
  FieldCell<Scal> fc;
  FieldCell<bool> fci;
  FieldCell<Vect> fcn;
  FieldCell<Scal> fck;
};

class Partstr : public TimerMesh {
 public:
  Partstr(Mesh& m) 
      : TimerMesh("partstr", m)
      , fc(m), fci(m, true), fcn(m), fca(m)
  {
    for (auto c : m.AllCells()) {
      Scal q = c.GetRaw();
      fc[c] = Rndn(q);
      fcn[c] = Vect(Rndn(q), Rndn(q), Rndn(q));
      fca[c] = Rndn(q);
    }
  }
  void F() override {
    volatile size_t ii = 0;

    ii = fcn[IdxCell(ii)][0];
  }

 private:
  FieldCell<Scal> fc;
  FieldCell<bool> fci;
  FieldCell<Vect> fcn;
  FieldCell<Scal> fca;
  FieldCell<Scal> fck;
};


// f=0: youngs
// f=1: height edim=3
// f=2: height edim=2
void Cmp(int f) {
  auto m = GetMesh(MIdx(8));
  FieldCell<Scal> fc(m);
  FieldCell<bool> fci(m, true);
  FieldCell<Vect> fcn(m);
  FieldCell<Vect> fcn2(m);
  FieldCell<Scal> fck(m);
  for (auto c : m.AllCells()) {
    fc[c] = Rnd(c.GetRaw());
  }
  if (f == 0) {
    Normal::CalcNormalYoung(m, fc, fci, fcn);
    Normal::CalcNormalYoung1(m, fc, fci, fcn2);
  } else if (f == 1) {
    size_t edim = 3;
    Normal::CalcNormalHeight(m, fc, fci, edim, true, fcn, fck);
    Normal::CalcNormalHeight1(m, fc, fci, edim, true, fcn2, fck);
  } else {
    size_t edim = 2;
    Normal::CalcNormalHeight(m, fc, fci, edim, true, fcn, fck);
    Normal::CalcNormalHeight1(m, fc, fci, edim, true, fcn2, fck);
  }

  Scal r = DiffMax(fcn, fcn2, m);
  Scal eps = 1e-15;
  if (r > eps) {
    std::vector<std::string> nn = 
        {"NormalYoung", "NormalHeight,edim=3", "NormalHeight,edim=2"};
    std::cerr
      << "Cmp " + nn[f] + ": max difference excceded "
      << std::scientific << std::setprecision(16)
      << r << " > " << eps << std::endl;
    std::terminate();
  }
}


// i: target index
// k: current index
// Output:
// create if k == i
// ++k
// p: pointer to new instance
template <class T>
void Try(Mesh& m, size_t i, size_t& k, Timer*& p) {
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
bool Run(const size_t i, Mesh& m,
         double& t, size_t& n, size_t& mem, std::string& name) {
  size_t k = 0;
  Timer* p = nullptr;

  Try<Young<0>>(m, i, k, p);
  Try<Young<1>>(m, i, k, p);
  Try<Height<0>>(m, i, k, p);
  Try<Height<1>>(m, i, k, p);
  Try<Partstr>(m, i, k, p);

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
  Cmp(0);
  Cmp(1);
  Cmp(2);

  // mesh size
  std::vector<MIdx> ss;
  for (int n : {8, 16, 32, 64}) {
    ss.emplace_back(n, n, 8);
  }

  const size_t ww = 16;

  using std::setw;

  for (auto s : ss) {
    auto m = GetMesh(s);
    const size_t nci = m.GetInBlockCells().size();
    std::cout 
        << "Mesh" 
        << " size=" << s
        << " incells=" << nci 
        << std::endl;

    int i = 0;
    double t;
    size_t n;
    size_t mem;
    std::string name;
    while (Run(i++, m, t, n, mem, name)) {
      std::cout 
          << setw(ww) << name 
          << setw(ww) << t * 1e9 / nci
          << setw(ww) << n
          << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }

}
