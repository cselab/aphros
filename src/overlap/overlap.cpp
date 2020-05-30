// Created by Petr Karnakov on 02.04.2019
// Copyright 2019 ETH Zurich

#include "overlap.hpp"

#include "overlap.h"

double GetSphereOverlap(
    const generic::Vect<double, 3>& x, const generic::Vect<double, 3>& h,
    const generic::Vect<double, 3>& c, double r) {
  std::array<vector_t, 8> vv{vector_t{-1, -1, -1}, vector_t{1, -1, -1},
                             vector_t{1, 1, -1},   vector_t{-1, 1, -1},
                             vector_t{-1, -1, 1},  vector_t{1, -1, 1},
                             vector_t{1, 1, 1},    vector_t{-1, 1, 1}};

  bool ai = true; // all inside sphere
  bool ao = true; // all outside sphere
  auto hh = h * 0.5;
  for (size_t i = 0; i < vv.size(); ++i) {
    auto& v = vv[i];
    v[0] = x[0] + v[0] * hh[0];
    v[1] = x[1] + v[1] * hh[1];
    v[2] = x[2] + v[2] * hh[2];
    bool in = generic::Vect<double, 3>(v[0], v[1], v[2]).sqrdist(c) < sqr(r);
    bool out = generic::Vect<double, 3>(v[0], v[1], v[2]).sqrdist(c) >
               sqr(r + h.max());
    ai = ai && in;
    ao = ao && out;
  }

  if (ai) {
    return 1.;
  }
  if (ao) {
    return 0.;
  }

  Hexahedron hex{vv};
  Sphere s{vector_t{c[0], c[1], c[2]}, r};

  return overlap(s, hex) / h.prod();
}
