#undef NDEBUG
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

// Echo Execute
#define EE(...); std::cerr << "\n" << #__VA_ARGS__ << std::endl; __VA_ARGS__;

Mesh GetMesh(MIdx s /*size in cells*/) {
  Rect<Vect> dom(Vect(0.1, 0.2, 0.1), Vect(1.1, 1.2, 1.3));
  MIdx b(-2, -3, -4); // lower index
  int hl = 2;         // halos 
  Vect doms = dom.GetDimensions();
  Vect h = dom.GetDimensions() / Vect(s);
  return InitUniformMesh<Mesh>(dom, b, s, hl, true, s);
}

class Empty : public Timer {
  Mesh& m;
 public:
  Empty(Mesh& m) : Timer("empty"), m(m) {}
  void F() override {}
};

class LoopPlain : public Timer {
  Mesh& m;
 public:
  LoopPlain(Mesh& m) : Timer("loop-plain"), m(m) {}
  void F() override {
    volatile size_t a = 0;
    for (size_t i = 0; i < m.GetNumCells(); ++i) {
      a += i;
    }
  }
};

class LoopAllCells : public Timer {
  Mesh& m;
 public:
  LoopAllCells(Mesh& m) : Timer("loop-allcells"), m(m) {}
  void F() override {
    volatile size_t a = 0;
    for (auto i : m.AllCells()) {
      a += i.GetRaw();
    }
  }
};

class LoopInCells : public Timer {
  Mesh& m;
 public:
  LoopInCells(Mesh& m) : Timer("loop-incells"), m(m) {}
  void F() override {
    volatile size_t a = 0;
    for (auto i : m.Cells()) {
      a += i.GetRaw();
    }
  }
};

class LoopAllFaces : public Timer {
  Mesh& m;
 public:
  LoopAllFaces(Mesh& m) : Timer("loop-allfaces"), m(m) {}
  void F() override {
    volatile size_t a = 0;
    for (auto i : m.AllFaces()) {
      a += i.GetRaw();
    }
  }
};

class LoopInFaces : public Timer {
  Mesh& m;
 public:
  LoopInFaces(Mesh& m) : Timer("loop-infaces"), m(m) {}
  void F() override {
    volatile size_t a = 0;
    for (auto i : m.Faces()) {
      a += i.GetRaw();
    }
  }
};


// Loop with array
class LoopArPlain : public Timer {
  Mesh& m;
  std::vector<Scal> v;
 public:
  LoopArPlain(Mesh& m) : Timer("loop-ar-plain"), m(m), v(m.GetNumCells()) {}
  void F() override {
    volatile Scal a = 0;
    for (size_t i = 0; i < m.GetNumCells(); ++i) {
      v[i] += a;
      a = v[i];
    }
  }
};

class LoopArAllCells : public Timer {
  Mesh& m;
  FieldCell<Scal> v;
 public:
  LoopArAllCells(Mesh& m) : Timer("loop-ar-allcells"), m(m), v(m) {}
  void F() override {
    volatile Scal a = 0;
    for (auto i : m.AllCells()) {
      v[i] += a;
      a = v[i];
    }
  }
};


class LoopMIdxAllCells : public Timer {
  Mesh& m;
 public:
  LoopMIdxAllCells(Mesh& m) : Timer("loop-midx-allcells"), m(m) {}
  void F() override {
    volatile size_t a = 0;
    for (auto i : m.AllCells()) {
      auto ii = m.GetBlockCells().GetMIdx(i);
      a += ii[0];
    }
  }
};

class LoopMIdxAllFaces : public Timer {
  Mesh& m;
 public:
  LoopMIdxAllFaces(Mesh& m) : Timer("loop-midx-allfaces"), m(m) {}
  void F() override {
    volatile size_t a = 0;
    for (auto i : m.AllFaces()) {
      auto ii = m.GetBlockFaces().GetMIdx(i);
      auto id = m.GetBlockFaces().GetDir(i);
      a += ii[0];
      a += id.GetLetter();
    }
  }
};

class Interp : public Timer {
  Mesh& m;
  FieldCell<Scal> fc;
  MapFace<std::shared_ptr<solver::CondFace>> mfc;
  FieldFace<Scal> ff;
 public:
  Interp(Mesh& m) : Timer("interp"), m(m), fc(m), ff(m) {
    for (auto i : m.AllCells()) {
      fc[i] = std::sin(i.GetRaw());
    }
    auto& bf = m.GetBlockFaces();
    for (auto i : m.Faces()) {
      if (bf.GetMIdx(i)[0] == 0 && bf.GetDir(i) == Dir::i) {
        mfc[i] = std::make_shared<solver::
            CondFaceGradFixed<Scal>>(0, 1);
      }
    }
    assert(mfc.size() > 0);
  }
  void F() override {
    volatile size_t a = 0;
    ff = solver::Interpolate(fc, mfc, m);
    a = ff[IdxFace(a)];
  }
};

class Grad : public Timer {
  Mesh& m;
  FieldCell<Vect> fc;
  FieldFace<Scal> ff;
 public:
  Grad(Mesh& m) : Timer("grad"), m(m), fc(m), ff(m) {
    for (auto i : m.AllFaces()) {
      ff[i] = std::sin(i.GetRaw());
    }
  }
  void F() override {
    static volatile size_t a = 0;
    fc = solver::Gradient(ff, m);
    a = fc[IdxCell(a % m.GetNumCells())][0];
  }
};


class ExplVisc : public Timer {
  Mesh& m;
  FieldCell<Vect> fcv;
  FieldCell<Vect> fcf;
  FieldFace<Scal> ffmu;
  MapFace<std::shared_ptr<solver::CondFace>> mfc;
  MapFace<std::shared_ptr<solver::CondFace>> mfcf;

 public:
  ExplVisc(Mesh& m) : Timer("explvisc"), m(m), fcv(m), fcf(m), ffmu(m) {
    for (auto i : m.AllCells()) {
      auto a = i.GetRaw();
      fcv[i] = Vect(std::sin(a), std::sin(a+1), std::sin(a+2));
    }
    for (auto i : m.AllFaces()) {
      auto a = i.GetRaw();
      ffmu[i] = std::sin(a);
    }
    auto& bf = m.GetBlockFaces();
    /*
    for (auto i : m.Faces()) {
      if (bf.GetMIdx(i)[0] == 0 && bf.GetDir(i) == Dir::i) {
        mfc[i] = std::make_shared<solver::
            CondFaceGradFixed<Scal>>(0);
      }
    }
    */
  }
  void F() override {
    static volatile size_t a = 0;
    using namespace solver;
    auto& mesh = m;
    for (size_t n = 0; n < dim; ++n) {
      FieldCell<Scal> fc = GetComponent(fcv, n);
      auto ff = Interpolate(fc, mfc, mesh);
      auto gc = Gradient(ff, mesh);
      auto gf = Interpolate(gc, mfcf, mesh); // adhoc: zero-der cond
      for (auto idxcell : mesh.SuCells()) {
        Vect sum = Vect::kZero;
        for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
          IdxFace idxface = mesh.GetNeighbourFace(idxcell, i);
          sum += gf[idxface] * (ffmu[idxface] * 
              mesh.GetOutwardSurface(idxcell, i)[n]);
        }
        fcf[idxcell] += sum / mesh.GetVolume(idxcell);
      }
    }
    a = fcf[IdxCell(a % m.GetNumCells())][0];
  }
};

// i: test index
// m: mesh
// Output:
// t: total per one call [sec] 
// n: number of calls 
// mem: memory usage in bytes
// name: test name
// Returns 1 if test with index i found
bool Run(size_t i, Mesh& m,
         double& t, size_t& n, size_t& mem, std::string& name) {
  size_t k = 0;
  Timer* p = nullptr;

  k++ == i && (p = new Empty(m));
  k++ == i && (p = new LoopPlain(m));
  k++ == i && (p = new LoopAllCells(m));
  k++ == i && (p = new LoopInCells(m));
  k++ == i && (p = new LoopAllFaces(m));
  k++ == i && (p = new LoopInFaces(m));
  k++ == i && (p = new LoopArPlain(m));
  k++ == i && (p = new LoopArAllCells(m));
  k++ == i && (p = new LoopMIdxAllCells(m));
  k++ == i && (p = new LoopMIdxAllFaces(m));
  k++ == i && (p = new Interp(m));
  k++ == i && (p = new Grad(m));
  k++ == i && (p = new ExplVisc(m));

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
  std::vector<MIdx> ss = {
        MIdx(4)
      , MIdx(8)
      , MIdx(16)
      , MIdx(32)
      , MIdx(64) 
      , MIdx(128) 
  };

  using std::setw;
  std::cout 
      << setw(20) << "name"
      << setw(20) << "t/allcells [ns]" 
      << setw(20) << "t/incells [ns]" 
      << setw(20) << "t [ns]" 
      << setw(20) << "iters"
      << setw(20) << "mem [MB]"
      << setw(20) << "mem/allcells [B]"
      << std::endl
      << std::endl;

  for (auto s : ss) {
    size_t mem0 = sysinfo::GetMem();
    auto m = GetMesh(s);
    const size_t nca = m.GetBlockCells().size();
    const size_t nci = m.GetInBlockCells().size();
    std::cout 
        << "Mesh" 
        << " size=" << s
        << " allcells=" << nca 
        << " incells=" << nci 
        << std::endl;

    int i = 0;
    double t;
    size_t n;
    size_t mem;
    std::string name;
    while (Run(i++, m, t, n, mem, name)) {
      size_t dmem = mem - mem0;
      std::cout 
          << setw(20) << name 
          << setw(20) << t * 1e9 / nca
          << setw(20) << t * 1e9 / nci
          << setw(20) << t * 1e9 * n
          << setw(20) << n
          << setw(20) << (dmem / double(1 << 20))
          << setw(20) << (dmem / nca)
          << std::endl;
    }
    std::cout << std::endl;
  }

}
