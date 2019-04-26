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
#include "timer.h"
#include "util/sysinfo.h"

const int dim = 3;
using MIdx = GMIdx<dim>;
using IdxCell = IdxCell;
using IdxFace = IdxFace;
using Dir = GDir<dim>;
using Scal = double;
using Vect = GVect<Scal, dim>;
using Mesh = MeshStructured<Scal, dim>;


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


class LoopPlain : public TimerMesh {
 public:
  LoopPlain(Mesh& m) : TimerMesh("loop-plain", m) {}
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
  LoopInCells(Mesh& m) : TimerMesh("loop-incells", m) {}
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
  LoopFldInCells(Mesh& m) : TimerMesh("loop-fld-incells", m), v(m) {}
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
  Interp(Mesh& m) : TimerMesh("interp", m), fc(m), ff(m) {
    for (auto i : m.AllCells()) {
      fc[i] = std::sin(i.GetRaw());
    }
    auto& bf = m.GetIndexFaces();
    for (auto i : m.Faces()) {
      if (bf.GetMIdx(i)[0] == 0 && bf.GetDir(i) == Dir::i) {
        mfc[i] = std::make_shared<solver::CondFaceGradFixed<Scal>>(0, 1);
      }
    }
    assert(mfc.size() > 0);
  }
  void F() override {
    volatile size_t a = 0;
    ff = solver::Interpolate(fc, mfc, m);
    a = ff[IdxFace(a)];
  }

 private:
  FieldCell<Scal> fc;
  MapFace<std::shared_ptr<solver::CondFace>> mfc;
  FieldFace<Scal> ff;
};

class Diffusion : public TimerMesh {
 public:
  Diffusion(Mesh& m) : TimerMesh("diffusion", m), fc(m), ff(m) {
    for (auto i : m.AllCells()) {
      fc[i] = std::sin(i.GetRaw());
    }
  }
  void F() override {
    volatile size_t ii = 0;

    // normal gradient
    for (auto f : m.Faces()) {
      auto cm = m.GetNeighbourCell(f, 0);
      auto cp = m.GetNeighbourCell(f, 1);
      ff[f] = (fc[cp] - fc[cm]);
    }

    // advance
    for (auto c : m.Cells()) {
      Scal s = 0;
      for (auto q : m.Nci(c)) {
        IdxFace f = m.GetNeighbourFace(c, q);
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
  DiffusionNeighb(Mesh& m) 
    : TimerMesh("diffusion-neighb", m), fc(m), ff(m) {
    for (auto i : m.AllCells()) {
      fc[i] = std::sin(i.GetRaw());
    }
  }
  void F() override {
    volatile size_t ii = 0;

    // advance
    for (auto c : m.Cells()) {
      Scal s = 0;
      for (auto q : m.Nci(c)) {
        IdxCell cc = m.GetNeighbourCell(c, q);
        s += fc[cc];
      }
      s += -6 * fc[c];
      fc[c] += s;
    }

    ii = fc[IdxCell(ii)];
  }

 private:
  FieldCell<Scal> fc;
  FieldFace<Scal> ff;
};

class DiffusionPlain : public TimerMesh {
 public:
  DiffusionPlain(Mesh& m) : TimerMesh("diffusion-plain", m), fc(m), fct(m) {
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
    const size_t hl = -wb[0];  // halo cells
    const size_t nx = ws[0];   // block size n^3
    const size_t ny = ws[1];   // block size n^3
    const size_t nz = ws[2];   // block size n^3

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
          size_t ixp = i + dx;
          size_t ixm = i - dx;
          size_t iyp = i + dy;
          size_t iym = i - dy;
          size_t izp = i + dz;
          size_t izm = i - dz;
          Scal s = (-6 * t[i] + 
              t[ixp] + t[ixm] + t[iyp] + t[iym] + t[izp] + t[izm]);
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

  Try<LoopPlain>(m, i, k, p);
  Try<LoopInCells>(m, i, k, p);
  //Try<LoopFldInCells>(m, i, k, p);
  //Try<Interp>(m, i, k, p);
  Try<Diffusion>(m, i, k, p);
  Try<DiffusionNeighb>(m, i, k, p);
  Try<DiffusionPlain>(m, i, k, p);

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
  for (int n : {32, 64, 128, 256, 512, 1024}) {
    ss.emplace_back(n, n, 8);
  }

  const size_t ww = 16;

  using std::setw;

  for (auto s : ss) {
    auto m = GetMesh(s);
    const size_t nca = m.GetAllBlockCells().size();
    const size_t nci = m.GetInBlockCells().size();
    std::cout 
        << "Mesh" 
        << " size=" << s
        << " allcells=" << nca 
        << " incells=" << nci 
        << std::endl;

    std::cout 
        << setw(ww) << "name"
        << setw(ww) << "t/allcells [ns]" 
        << setw(ww) << "t/incells [ns]" 
        << setw(ww) << "t [ns]" 
        << setw(ww) << "iters"
        << std::endl;
    for (size_t q = 0; q < ww * 5; ++q) {
      std::cout << "-";
    }
    std::cout << std::endl;

    int i = 0;
    double t;
    size_t n;
    size_t mem;
    std::string name;
    while (Run(i++, m, t, n, mem, name)) {
      std::cout 
          << setw(ww) << name 
          << setw(ww) << t * 1e9 / nca
          << setw(ww) << t * 1e9 / nci
          << setw(ww) << t * 1e9 * n
          << setw(ww) << n
          << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }

}
