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
      : TimerMesh("young-orig" + std::to_string(id), m)
      , fc(m), fci(m, true)
  {
    for (auto i : m.AllCells()) {
      fc[i] = std::sin(i.GetRaw()*0.17);
    }
  }
  void F() override {
    volatile size_t ii = 0;

    if (id == 0) {
      Normal::CalcNormalYoung(m, fc, fci, fcn);
    } else {
      Normal::CalcNormalYoung2(m, fc, fci, fcn);
    }

    ii = fcn[IdxCell(ii)][0];
  }

 private:
  FieldCell<Scal> fc;
  FieldCell<bool> fci;
  FieldCell<Vect> fcn;
};


class Grad : public TimerMesh {
 public:
  Grad(Mesh& m) 
      : TimerMesh("grad", m)
      , fc(m)
  {
    for (auto i : m.AllCells()) {
      fc[i] = std::sin(i.GetRaw());
    }
  }
  void F() override {
    volatile size_t ii = 0;

    fcn.Reinit(m);

    Normal::CalcNormalGrad(m, fc, 
        MapFace<std::shared_ptr<solver::CondFace>>(), fcn);

    ii = fcn[IdxCell(ii)][0];
  }

 private:
  FieldCell<Scal> fc;
  FieldCell<Vect> fcn;
};

void Cmp() {
    auto m = GetMesh(MIdx(8));
    FieldCell<Scal> fc(m);
    FieldCell<bool> fci(m, true);
    FieldCell<Vect> fcn(m);
    for (auto c : m.AllCells()) {
      fc[c] = std::sin(std::sin(c.GetRaw()) * 123);
    }
    for (size_t q : {0, 1}) {
      if (q == 0) Normal::CalcNormalYoung(m, fc, fci, fcn);
      if (q == 1) Normal::CalcNormalYoung2(m, fc, fci, fcn);
      size_t i = 3;
      for (auto c : m.Cells()) {
        std::cout << fcn[c] << std::endl;
        if (!--i) break;
      }
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
  Try<Grad>(m, i, k, p);

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
