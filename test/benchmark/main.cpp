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
    for (size_t i = 0; i < m.GetNumCells(); ++i) {
      volatile size_t ii = i;
    }
  }
};

int main() {
  std::vector<MIdx> ss = {
      MIdx(8, 8, 8), MIdx(32, 32, 32),
      MIdx(8, 8, 2), MIdx(32, 32, 2),
      MIdx(8, 8, 1), MIdx(32, 32, 1)
    };

  using std::setw;
  std::cout 
      << setw(15) 
      << "name"
      << setw(15) 
      << "total [ns]" 
      << setw(15) 
      << "percell [ns]" 
      << setw(15) 
      << "iters"
      << std::endl
      << std::endl;

  for (auto s : ss) {
    std::cout << "Mesh " << s << "" << std::endl;
    auto m = GetMesh(s);

    std::vector<Timer*> tt {
          new Empty(m)
        , new LoopPlain(m)
      };

    for (auto& t : tt) {
      auto e = t->Run();
      std::cout 
          << setw(15) 
          << t->GetName()
          << setw(15) 
          << e.first * 1e9 
          << setw(15) 
          << e.first * 1e9 / m.GetNumCells() 
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
