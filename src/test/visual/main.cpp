// Created by Petr Karnakov on 08.03.2021
// Copyright 2021 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>

#include "geom/mesh.h"
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
  return InitUniformMesh<M>(
      domain, MIdx(0), size, halos, true, true, size, 0);
}

void TestRender() {
  using U = util::Visual<M>;
  typename U::Canvas canvas(MIdx(64));
  auto m = GetMesh(MIdx(64));
  typename U::CanvasView view(canvas);
  using Float3 = typename U::Float3;
  FieldCell<Float3> fc_color(m, Float3(0));
  FieldCell<Scal> fc(m, 0);
  for (auto c : m.CellsM()) {
    fc_color[c][0] = Vect(0, 0).dist(c.center);
    fc_color[c][1] = Vect(0, 0.5).dist(c.center);
    fc_color[c][2] = Vect(0.5, 0).dist(c.center);
    fc[c] = Vect(0.5, 0.5).dist(c.center);
  }
  typename U::Colormap cmap;
  cmap.min = 0;
  cmap.max = 1;
  cmap.color_min = Float3(1, 0, 0);
  cmap.color_max = Float3(0, 1, 0);
  cmap.opacity_min = 1;
  cmap.opacity_max = 0;
  U::AppendField(fc_color, fc, cmap, m);

  U::RenderColorField(fc_color, view, m);
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
