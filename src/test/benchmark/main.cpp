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
  return InitUniformMesh<Mesh>(dom, b, s, hl, true, s);
}

class TimerMesh : public Timer {
 public:
  TimerMesh(const std::string& name, Mesh& m) : Timer(name), m(m) {}

 protected:
  Mesh& m;
};

/*
class IFactoryTimerMesh {
 public:
  virtual TimerMesh* Make(Mesh& m) = 0;
};

std::vector<std::unique_ptr<IFactoryTimerMesh>> reg;

template <class T>
class FactoryTimerMesh : IFactoryTimerMesh {
 public:
  TimerMesh* Make(Mesh& m) override {
    return new T(m);
  }
};
*/


class Empty : public TimerMesh {
 public: 
  Empty(Mesh& m) : TimerMesh("empty", m) {}
  void F() override {}
};

class LoopPlain : public TimerMesh {
 public:
  LoopPlain(Mesh& m) : TimerMesh("loop-plain", m) {}
  void F() override {
    volatile size_t a = 0;
    for (size_t i = 0; i < m.GetAllBlockCells().GetSize().prod(); ++i) {
      a += i;
    }
  }
};

class LoopAllCells : public TimerMesh {
 public:
  LoopAllCells(Mesh& m) : TimerMesh("loop-allcells", m) {}
  void F() override {
    volatile size_t a = 0;
    for (auto i : m.AllCells()) {
      a += i.GetRaw();
    }
  }
};

class LoopInCells : public TimerMesh {
 public:
  LoopInCells(Mesh& m) : TimerMesh("loop-incells", m) {}
  void F() override {
    volatile size_t a = 0;
    for (auto i : m.Cells()) {
      a += i.GetRaw();
    }
  }
};

class LoopAllFaces : public TimerMesh {
 public:
  LoopAllFaces(Mesh& m) : TimerMesh("loop-allfaces", m) {}
  void F() override {
    volatile size_t a = 0;
    for (auto i : m.AllFaces()) {
      a += i.GetRaw();
    }
  }
};

class LoopInFaces : public TimerMesh {
 public:
  LoopInFaces(Mesh& m) : TimerMesh("loop-infaces", m) {}
  void F() override {
    volatile size_t a = 0;
    for (auto i : m.Faces()) {
      a += i.GetRaw();
    }
  }
};


// Loop access to field
class LoopFldPlain : public TimerMesh {
 public:
  LoopFldPlain(Mesh& m) : TimerMesh("loop-fld-plain", m), v(GRange<IdxCell>(m).size()) {}
  void F() override {
    volatile Scal a = 0;
    for (size_t i = 0; i < v.size(); ++i) {
      v[i] += a;
      a = v[i];
    }
  }

 private:
  std::vector<Scal> v;
};

class LoopFldAllCells : public TimerMesh {
 public:
  LoopFldAllCells(Mesh& m) : TimerMesh("loop-fld-allcells", m), v(m) {}
  void F() override {
    volatile Scal a = 0;
    for (auto i : m.AllCells()) {
      v[i] += a;
      a = v[i];
    }
  }

 private:
  FieldCell<Scal> v;
};


class LoopMIdxAllCells : public TimerMesh {
 public:
  LoopMIdxAllCells(Mesh& m) : TimerMesh("loop-midx-allcells", m) {}
  void F() override {
    volatile size_t a = 0;
    for (auto c : m.AllCells()) {
      auto w = m.GetIndexCells().GetMIdx(c);
      a += w[0];
    }
  }
};

class LoopMIdxAllFaces : public TimerMesh {
 public:
  LoopMIdxAllFaces(Mesh& m) : TimerMesh("loop-midx-allfaces", m) {}
  void F() override {
    volatile size_t a = 0;
    for (auto f : m.AllFaces()) {
      auto wd = m.GetIndexCells().GetMIdxDir(f);
      a += wd.first[0];
      a += wd.second.GetLetter();
    }
  }
};

class CellVolume : public TimerMesh {
 public:
  CellVolume(Mesh& m) : TimerMesh("m-cell-volume", m) {}
  void F() override {
    volatile size_t a = 1;
    for (auto c : m.Cells()) {
      a += m.GetVolume(c);
    }
  }
};

class CellCenter : public TimerMesh {
 public:
  CellCenter(Mesh& m) : TimerMesh("m-cell-center", m) {}
  void F() override {
    volatile int a = 0;
    for (auto c : m.Cells()) {
      a += m.GetCenter(c)[0];
    }
  }
};

class CellNCell : public TimerMesh {
 public:
  CellNCell(Mesh& m) : TimerMesh("m-cell-n-cell", m) {}
  void F() override {
    volatile size_t a = 1;
    for (auto c : m.Cells()) {
      for (auto q : m.Nci(c)) {
        a += m.GetNeighbourCell(c, q).GetRaw();
      }
    }
  }
};

class CellNFace : public TimerMesh {
 public:
  CellNFace(Mesh& m) : TimerMesh("m-cell-n-face", m) {}
  void F() override {
    volatile size_t a = 1;
    for (auto c : m.Cells()) {
      for (auto q : m.Nci(c)) {
        a += m.GetNeighbourFace(c, q).GetRaw();
      }
    }
  }
};

class CellOutward : public TimerMesh {
 public:
  CellOutward(Mesh& m) : TimerMesh("m-cell-outward", m) {}
  void F() override {
    volatile size_t a = 0;
    for (auto c : m.Cells()) {
      for (auto q : m.Nci(c)) {
        a += m.GetOutwardFactor(c, q);
      }
    }
  }
};


class CellNNode : public TimerMesh {
 public:
  CellNNode(Mesh& m) : TimerMesh("m-cell-n-node", m) {}
  void F() override {
    volatile size_t a = 1;
    for (auto c : m.Cells()) {
      for (size_t q = 0; q < m.GetNumNeighbourNodes(c); ++q) {
        a += m.GetNeighbourNode(c, q).GetRaw();
      }
    }
  }
};

class FaceCenter : public TimerMesh {
 public:
  FaceCenter(Mesh& m) : TimerMesh("m-face-center", m) {}
  void F() override {
    volatile int a = 0;
    for (auto f : m.Faces()) {
      a += m.GetCenter(f)[0];
    }
  }
};

class FaceSurf : public TimerMesh {
 public:
  FaceSurf(Mesh& m) : TimerMesh("m-face-surf", m) {}
  void F() override {
    volatile int a = 0;
    for (auto f : m.Faces()) {
      a += m.GetSurface(f)[0];
    }
  }
};

class FaceArea : public TimerMesh {
 public:
  FaceArea(Mesh& m) : TimerMesh("m-face-area", m) {}
  void F() override {
    volatile int a = 0;
    for (auto f : m.Faces()) {
      a += m.GetSurface(f)[0];
    }
  }
};

class FaceNCell : public TimerMesh {
 public:
  FaceNCell(Mesh& m) : TimerMesh("m-face-n-cell", m) {}
  void F() override {
    volatile size_t a = 1;
    for (auto f : m.Faces()) {
      for (size_t q = 0; q < m.GetNumNeighbourCells(f); ++q) {
        a += m.GetNeighbourCell(f, q).GetRaw();
      }
    }
  }
};

class FaceNNode : public TimerMesh {
 public:
  FaceNNode(Mesh& m) : TimerMesh("m-face-n-node", m) {}
  void F() override {
    volatile size_t a = 1;
    for (auto f : m.Faces()) {
      for (size_t q = 0; q < m.GetNumNeighbourNodes(f); ++q) {
        a += m.GetNeighbourNode(f, q).GetRaw();
      }
    }
  }
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

 private:
  FieldCell<Scal> fc;
  MapFace<std::shared_ptr<solver::CondFace>> mfc;
  FieldFace<Scal> ff;
};

class Grad : public TimerMesh {
 public:
  Grad(Mesh& m) : TimerMesh("grad", m), fc(m), ff(m) {
    for (auto i : m.AllFaces()) {
      ff[i] = std::sin(i.GetRaw());
    }
  }
  void F() override {
    static volatile size_t a = 0;
    fc = solver::Gradient(ff, m);
    a += fc[IdxCell(0)][0];
  }

 private:
  FieldCell<Vect> fc;
  FieldFace<Scal> ff;
};


class ExplVisc : public TimerMesh {
 public:
  ExplVisc(Mesh& m) : TimerMesh("explvisc", m), fcv(m), fcf(m), ffmu(m) {
    for (auto i : m.AllCells()) {
      auto a = i.GetRaw();
      fcv[i] = Vect(std::sin(a), std::sin(a+1), std::sin(a+2));
    }
    for (auto i : m.AllFaces()) {
      auto a = i.GetRaw();
      ffmu[i] = std::sin(a);
    }
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
    a += fcf[IdxCell(0)][0];
  }

 private:
  FieldCell<Vect> fcv;
  FieldCell<Vect> fcf;
  FieldFace<Scal> ffmu;
  MapFace<std::shared_ptr<solver::CondFace>> mfc;
  MapFace<std::shared_ptr<solver::CondFace>> mfcf;
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
    const size_t nca = m.GetAllBlockCells().size();
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
