// Created by Petr Karnakov on 08.03.2021
// Copyright 2021 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>

#include "geom/mesh.h"
#include "parse/codeblocks.h"
#include "parse/parser.h"
#include "parse/vars.h"
#include "util/visual.h"

const int dim = 2;
using M = MeshStructured<double, dim>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;
using Vect3 = generic::Vect<Scal, 3>;

M GetMesh(MIdx size) {
  const Rect<Vect> domain(Vect(0), Vect(1));
  const size_t halos = 1;
  return InitUniformMesh<M>(domain, MIdx(0), size, halos, true, true, size, 0);
}

void TestRender() {
  using U = util::Visual<M>;
  typename U::Canvas canvas(MIdx(256));
  auto m = GetMesh(MIdx(64));
  typename U::CanvasView view(canvas);
  using Float3 = typename U::Float3;
  FieldCell<Float3> fc_color(m, Float3(0));
  FieldCell<Scal> fc(m, 0);
  FieldCell<Scal> fc2(m, 0);
  for (auto c : m.CellsM()) {
    fc[c] = Vect(0.5, 0.5).dist(c.center);
    fc2[c] = Vect(0.2, 0.8).dist(c.center);
  }

  auto entries = U::ParseEntries("vislist");

  auto get_field = [&](std::string name) -> FieldCell<Scal> {
    if (name == "p") {
      return fc;
    }
    if (name == "vx") {
      return fc2;
    }
    fassert(false, "Unknown field '" + name + "'");
  };

  U::RenderEntriesToField(fc_color, entries, get_field, m);

  U::RenderToCanvasNearest(view, fc_color, m);
  const auto path = "out.ppm";
  const auto path2 = "out2.ppm";
  U::WritePpm(path, view);
  auto canvas2 = U::ReadPpm(path);
  auto view2 = U::CanvasView(canvas2);
  U::WritePpm(path2, view2);
}

int main() {
  TestRender();
}
