// Created by Petr Karnakov on 09.07.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <vector>

#include "geom/vect.h"

template <class Scal>
struct UHeight {
  // Computes the height if well-defined in a cell close to the interface.
  // V: vector type with size() and operator[]
  // uu: column of volume fractions
  // n: component of normal in the direction of column
  // Returns:
  // offset from cell center to the interface, and NaN if not well-defined.
  template <class V>
  static Scal Good(const V& uu, Scal n) {
    auto I = [](Scal a) { return a > 0 && a < 1; }; // true if interface

    const Scal nan = GetNan<Scal>();

    const size_t si = uu.size();
    const size_t sih = si / 2;
    size_t icenter = sih; // closest interface to center
    while (icenter < si && !I(uu[icenter])) {
      if (icenter > sih) {
        icenter = si - icenter - 1;
      } else {
        icenter = si - icenter;
      }
    }
    if (icenter >= si) {
      return nan;
    }

    size_t im = icenter; // closest pure cell below
    while (im < si && I(uu[im])) {
      --im;
    }
    if (im >= si) {
      return nan;
    }

    size_t ip = icenter; // closest pure cell above
    while (ip < si && I(uu[ip])) {
      ++ip;
    }
    if (ip >= si) {
      return nan;
    }

    if ((uu[ip] - uu[im]) * n < 0) {
      Scal u = 0;
      for (size_t i = im; i <= ip; ++i) {
        u += uu[i];
      }
      return u + im * uu[im] + (si - ip - 1) * uu[ip];
    }
    return nan;
  }

  // Computes the height if well-defined in an interfacial cell.
  // V: vector type with size() and operator[]
  // uu: column of volume fractions
  // Returns:
  // offset from cell center to the interface, and NaN if not well-defined.
  template <class V>
  static Scal Good(const V& uu) {
    auto I = [](Scal a) { return a > 0 && a < 1; }; // true if interface

    const Scal nan = GetNan<Scal>();

    const size_t si = uu.size();
    const size_t sih = si / 2;

    const size_t icenter = sih; // center
    if (!I(uu[icenter])) { // check center is interface
      return nan;
    }

    size_t im = icenter; // closest pure cell below
    while (im < si && I(uu[im])) {
      --im;
    }
    if (im >= si) {
      return nan;
    }

    size_t ip = icenter; // closest pure cell above
    while (ip < si && I(uu[ip])) {
      ++ip;
    }
    if (ip >= si) {
      return nan;
    }

    if (uu[ip] != uu[im]) {
      Scal u = 0.; // sum over column
      for (size_t i = im; i <= ip; ++i) {
        u += uu[i];
      }
      // add missing pure cells
      u += im * uu[im] + (si - ip - 1) * uu[ip];

      // assume cell center is x=0,
      // if uu[im] == 1, u is distance from x=-si*0.5
      // if uu[im] == 0, u is distance from x=+si*0.5

      if (uu[im] == 0) {
        u = si - u;
      }

      // u is distance from x=-(sih+0.5)

      // subtract center
      u -= si * 0.5;

      // u is distance from x=0
      return u;
    }
    return nan;
  }
};
