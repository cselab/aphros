// Created by Petr Karnakov on 03.11.2020
// Copyright 2020 ETH Zurich

#pragma once
#include <x86intrin.h>

namespace util {

std::ostream& operator<<(std::ostream& out, const __m256d& d) {
  constexpr int width = sizeof(__m256d) / sizeof(double);
  double dd[width];
  _mm256_storeu_pd(dd, d);
  bool first = true;
  for (size_t i = 0; i < width; ++i) {
    if (!first) {
      out << " ";
    }
    first = false;
    out << dd[i];
  }
  return out;
}

struct Soa {
  static void ToAos(
      const __m256d& x, const __m256d& y, const __m256d& z, __m256d& v0,
      __m256d& v1, __m256d& v2) {
    // x: x0 x1 x2 x3
    // y: y0 y1 y2 y3
    // z: z0 z1 z2 z3
    auto xp = _mm256_permute2f128_pd(x, x, 0b00000001); // x2 x3 x0 x1
    auto zp = _mm256_permute2f128_pd(z, z, 0b00000001); // z2 z3 z0 z1
    auto xs = _mm256_blend_pd(x, xp, 0b1010); /// x0 x3 x2 x1
    auto ys = _mm256_shuffle_pd(y, y, 0b0101); // y1 y0 y3 y2
    auto zs = _mm256_blend_pd(z, zp, 0b0101); /// z2 z1 z0 z3
    auto xyb = _mm256_blend_pd(xs, ys, 0b1010); // x0 y0 x2 y2
    auto yzb = _mm256_blend_pd(ys, zs, 0b1010); // y1 z1 y3 z3
    auto zxb = _mm256_blend_pd(zs, xs, 0b1010); // z2 x3 z0 x1
    v0 = _mm256_blend_pd(xyb, zxb, 0b1100); // x0 y0 z0 x1
    v1 = _mm256_blend_pd(yzb, xyb, 0b1100); // y1 z1 x2 y2
    v2 = _mm256_blend_pd(zxb, yzb, 0b1100); // y2 z3 x3 y3
  }

  static void StoreAsAos(
      const __m256d& x, const __m256d& y, const __m256d& z, double* mem) {
    __m256d d0;
    __m256d d1;
    __m256d d2;

    ToAos(x, y, z, d0, d1, d2);

    _mm256_storeu_pd(mem + 0, d0);
    _mm256_storeu_pd(mem + 4, d1);
    _mm256_storeu_pd(mem + 8, d2);
  };

  static void Normalize1(__m256d& x, __m256d& y, __m256d& z) {
    auto xa = _mm256_max_pd(_mm256_sub_pd(_mm256_setzero_pd(), x), x);
    auto ya = _mm256_max_pd(_mm256_sub_pd(_mm256_setzero_pd(), y), y);
    auto za = _mm256_max_pd(_mm256_sub_pd(_mm256_setzero_pd(), z), z);
    auto sum = _mm256_add_pd(xa, ya);
    sum = _mm256_add_pd(sum, za);
    sum = _mm256_sub_pd(_mm256_setzero_pd(), sum);
    x = _mm256_div_pd(x, sum);
    y = _mm256_div_pd(y, sum);
    z = _mm256_div_pd(z, sum);
  };
};

} // namespace util
