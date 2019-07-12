#pragma once

#include <vector>

#include "geom/vect.h"

template <class Scal>
struct UHeight {
  // checks that the heiht in column is well-defined
  // V: vector type with size() and operator[]
  // uu: column of volume fractions
  // n: component of normal in the direction of column
  template <class V>
  static Scal Good(const V& uu, Scal n) {
    auto I = [](Scal a) { return a > 0 && a < 1; }; // true if interface

    const Scal nan = GetNan<Scal>();

    const size_t si = uu.size();
    const size_t sih = si / 2;
    size_t i = sih; // closest interface to center
    while (i < si && !I(uu[i])) {
      if (i > sih) {
        i = si - i - 1;
      } else {
        i = si - i;
      }
    }
    if (i >= si) { return nan; }

    size_t im = i; // closest pure cell below
    while (im < si && I(uu[im])) { --im; }
    if (im >= si) { return nan; }

    size_t ip = i; // closest pure cell above
    while (ip < si && I(uu[ip])) { ++ip; }
    if (ip >= si) { return nan; }

    if ((uu[ip] - uu[im]) * n < 0) {
      Scal u = 0;
      for (size_t i = im; i <= ip; ++i) {
        u += uu[i];
      }
      return u + im * uu[im] + (si - ip - 1) * uu[ip];
    }
    return nan;
  }

  // checks that the heiht in column is well-defined
  // V: vector type with size() and operator[]
  // uu: column of volume fractions
  // Returns:
  // offset from cell center to the interface
  template <class V>
  static Scal Good(const V& uu) {
    auto I = [](Scal a) { return a > 0 && a < 1; }; // true if interface

    const Scal nan = GetNan<Scal>();

    const size_t si = uu.size();
    const size_t sih = si / 2;

    size_t i = sih; // center
    if (!I(uu[i])) { // check center is interface
      return nan;
    }

    size_t im = i; // closest pure cell below
    while (im < si && I(uu[im])) { --im; }
    if (im >= si) { return nan; }

    size_t ip = i; // closest pure cell above
    while (ip < si && I(uu[ip])) { ++ip; }
    if (ip >= si) { return nan; }

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
