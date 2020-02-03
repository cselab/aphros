// Created by Petr Karnakov on 03.02.2020
// Copyright 2020 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>

#include "geom/mesh.h"
#include "solver/approx.h"
#include "solver/cond.h"
#include "solver/solver.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

// h: cell size
M GetMesh(Scal h) {
  const MIdx size(32);
  Rect<Vect> dom(Vect(0), Vect(h) * Vect(size));
  return InitUniformMesh<M>(dom, MIdx(0), size, 2, true, true, size, 0);
}

int main() {
}
