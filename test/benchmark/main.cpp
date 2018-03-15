#include <sstream>
#include <iostream>
#include <cassert>
#include <functional>
#include <cmath>
#include <string>
#include <memory>
#include "timer.h"

#include "mesh.hpp"
#include "mesh3d.hpp"
#include "solver.hpp"
#include "metrics.hpp"

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
  Empty(Mesh& m)
    : Timer("empty"), m(m)
  {}
  double F() override {
    int a = 0;
    for (auto i : m.Cells()) {
      a += 1.;
    }
    return a;
  }
};

int main() {
  std::vector<MIdx> ss = {
      MIdx(8, 8, 8), MIdx(16, 16, 16),
      MIdx(8, 8, 2), MIdx(16, 16, 2),
      MIdx(8, 8, 1), MIdx(16, 16, 1)
    };

  for (auto s : ss) {
    std::cout << "*** Mesh " << s << " ***" << std::endl;
    auto m = GetMesh(s);

    std::vector<Timer*> tt {
        new Empty(m)
      };

    for (auto& t : tt) {
      auto e = t->Run();
      std::cout 
          << t->GetName() << "\t "
          << "total=" << e.first * 1e6 << " us \t"
          << "per-cell=" << e.first * 1e6 / m.GetNumCells() << " us \t"
          << "iter=" << e.second << " "
          << std::endl;
    }

    for (auto& t : tt) {
      delete t;
    }
    std::cout << std::endl;
  }

}
