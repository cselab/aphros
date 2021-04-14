// Created by Petr Karnakov on 07.10.2018
// Copyright 2018 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "debug/isnan.h"
#include "geom/field.h"
#include "geom/mesh.h"
#include "geom/range.h"
#include "solver/sphavg.h"
#include "solver/trackerm.h"

// try-catch
#define TR(...)                           \
  try {                                   \
    __VA_ARGS__                           \
  } catch (const std::runtime_error& e) { \
    std::cout << e.what() << std::endl;   \
  }

// print current func
#define PF std::cout << std::endl << __func__ << std::endl;

const int dim = 3;
using MIdx = generic::MIdx<dim>;
using Dir = GDir<dim>;
using Scal = double;
using Vect = generic::Vect<Scal, dim>;
using M = MeshCartesian<Scal, dim>;

M GetMesh() {
  Rect<Vect> dom(Vect(0., 0., 0.), Vect(1., 1., 1.));
  MIdx b(0, 0, 0); // lower index
  MIdx s(10, 10, 10); // size in cells
  int hl = 2; // halos
  return {b, s, dom, hl, true, true, MIdx(20, 20, 20), 0};
}

void TestAvg() {
  PF;

  M m = GetMesh();
  Sphavg<M> sa(m, 2);

  FieldCell<Scal> fcu(m);
  FieldCell<Vect> fcv(m);
  FieldCell<Vect> fcvm(m);
  FieldCell<Scal> fcp(m);
  Scal dt = 1.;
  using Sph = typename Sphavg<M>::Sph;
  std::vector<Sph> ss;
  ss.emplace_back(Vect(1.9), 0.2, 0.05);

  sa.Update(fcu, fcv, fcvm, dt, fcp, ss);
  sa.Update(fcu, fcv, fcvm, dt, fcp, ss);
  sa.Update(fcu, fcv, fcvm, dt, fcp, ss);
  sa.Update(fcu, fcv, fcvm, dt, fcp, ss);
}

void TestPack() {
  using R = Trackerm<M>;
  auto p = [](int w0, int w1, int w2) {
    MIdx w(w0, w1, w2);
    std::cout << w << " " << R::Unpack(R::Pack(w)) << std::endl;
  };
  std::cout << "sizeof(Bit)=" << sizeof(typename R::Bit) << std::endl;
  p(1, 2, 3);
  p(30000, 31000, 32000);
  p(-30000, -31000, -32000);
  p(60000, 62000, 64000);
}

int main() {
  TestAvg();
  TestPack();
}
