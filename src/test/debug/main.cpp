// Created by Petr Karnakov on 03.10.2018
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
#include "util/timer.h"

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
using MIdx = GMIdx<dim>;
using Dir = GDir<dim>;
using Scal = double;
using Vect = generic::Vect<Scal, dim>;
using M = MeshStructured<Scal, dim>;

M GetMesh() {
  Rect<Vect> dom(Vect(0., 0., 0.), Vect(1., 1., 1.));
  MIdx b(0, 0, 0); // lower index
  MIdx s(2, 2, 1); // size in cells
  int hl = 0; // halos
  return InitUniformMesh<M>(dom, b, s, hl, true, true, s, 0);
}

void TestRange() {
  PF;

  GRange<size_t> r(0, 10);
  GField<Scal, size_t> f(r);
  f[0] = GetNan<Scal>();
  TR(CheckNan(f, "f", r);)
}

void TestMesh() {
  PF;

  M m = GetMesh();
  m.flags.check_nan = true;
  FieldCell<Scal> f(m);
  f[IdxCell(0)] = GetNan<Scal>();
  TR(CheckNan(f, "f", m);)
  TR(CHECKNAN(f, m.flags.check_nan);)
}

void TestAssert() {
  TR(fassert(1);)
  TR(fassert(1, "message");)
  TR(fassert(0);)
  TR(fassert(0, "message");)
  TR(fassert_equal(0, 0);)
  TR(fassert_equal(0, 1);)
  const std::string a = "a";
  const std::string b = "b";
  TR(fassert_equal(a, a);)
  TR(fassert_equal(a, b);)
}

int main() {
  TestRange();
  TestMesh();
  TestAssert();
}
