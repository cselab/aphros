#include <sstream>
#include <iostream>
#include <cassert>
#include <functional>
#include <cmath>
#include <string>
#include <memory>
#include "timer.h"
#include <iomanip>

#include "mesh.hpp"
#include "mesh3d.hpp"
#include "solver.hpp"

const int dim = 3;
using MIdx = geom::GMIdx<dim>;
using IdxCell = geom::IdxCell;
using IdxFace = geom::IdxFace;
using Dir = geom::GDir<dim>;
using Scal = double;
using Vect = geom::GVect<Scal, dim>;
using Mesh = geom::MeshStructured<Scal, dim>;

// Echo Execute
#define EE(...); std::cerr << "\n" << #__VA_ARGS__ << std::endl; __VA_ARGS__;

Mesh GetMesh(MIdx s /*size in cells*/) {
  geom::Rect<Vect> dom(Vect(0.1, 0.2, 0.1), Vect(1.1, 1.2, 1.3));
  MIdx b(-2, -3, -4); // lower index
  int hl = 2;         // halos 
  Vect doms = dom.GetDimensions();
  Vect h = dom.GetDimensions() / Vect(s);
  return geom::InitUniformMesh<Mesh>(dom, b, s, hl);
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
  geom::FieldCell<Scal> v;
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
  geom::FieldCell<Scal> fc;
  geom::MapFace<std::shared_ptr<solver::ConditionFace>> mfc;
  geom::FieldFace<Scal> ff;
 public:
  Interp(Mesh& m) : Timer("interp"), m(m), fc(m), ff(m) {
    for (auto i : m.AllCells()) {
      fc[i] = std::sin(i.GetRaw());
    }
    auto& bf = m.GetBlockFaces();
    for (auto i : m.Faces()) {
      if (bf.GetMIdx(i)[0] == 0 && bf.GetDir(i) == Dir::i) {
        mfc[i] = std::make_shared<solver::
            ConditionFaceDerivativeFixed<Scal>>(0);
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
  geom::FieldCell<Vect> fc;
  geom::FieldFace<Scal> ff;
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



int main() {
  std::vector<MIdx> ss = {
      //MIdx(16, 16, 16), MIdx(32, 32, 32), MIdx(64, 64, 64),
      //MIdx(16, 16, 1), MIdx(32, 32, 1), MIdx(64, 64, 1),
        MIdx(4)
      , MIdx(8)
      , MIdx(16)
      , MIdx(32)
      , MIdx(64) 
  };

  using std::setw;
  std::cout 
      << setw(20) 
      << "name"
      << setw(15) 
      << "percell [ns]" 
      << setw(15) 
      << "total [ns]" 
      << setw(15) 
      << "iters"
      << std::endl
      << std::endl;

  for (auto s : ss) {
    auto m = GetMesh(s);
    std::cout 
        << "Mesh" 
        << " size=" << s
        << " numcells=" << m.GetNumCells() 
        << std::endl;

    std::vector<Timer*> tt {
          new Empty(m)
        , new LoopPlain(m)
        , new LoopAllCells(m)
        , new LoopInCells(m)
        , new LoopAllFaces(m)
        , new LoopInFaces(m)
        , new LoopArPlain(m)
        , new LoopArAllCells(m)
        , new LoopMIdxAllCells(m)
        , new LoopMIdxAllFaces(m)
        , new Interp(m)
        , new Grad(m)
      };

    for (auto& t : tt) {
      auto e = t->Run();
      std::cout 
          << setw(20) 
          << t->GetName()
          << setw(15) 
          << e.first * 1e9 / m.GetNumCells() 
          << setw(15) 
          << e.first * 1e9 
          << setw(15) 
          << e.second
          << std::endl;
    }

    for (auto& t : tt) {
      delete t;
    }
    std::cout << std::endl;
  }

}
