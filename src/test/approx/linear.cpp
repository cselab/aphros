// Created by Petr Karnakov on 07.02.2020
// Copyright 2020 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <random>
#include <sstream>
#include <array>

#include "solver/embed.h"
#include "geom/vect.h"

using Scal = double;
using Vect = generic::Vect<Scal, 3>;

void TestSolveLinear() {
  constexpr size_t N = 4;
  std::array<Scal, N * N> a = {
    0, 0, 0.1, 9.8,
    0.7, 1, 0.8, 0.3,
    1, 0, 0.3, 0.2,
    1, 0, 0.3, 0.1,
  };
  std::array<Scal, N> b = {
    1, 2, 8,
  };
  auto x = SolveLinear(a, b);
  using V = generic::Vect<Scal, N>;
  std::cout << "x=" << V(x) << std::endl;
  std::cout << "b=" << V(b) << std::endl;
  std::cout << "a*x=" << V(Mul(a, x)) << std::endl;
}

void TestFitLinear(const std::vector<Vect>& xx, const std::vector<Scal>& uu) {
  std::cout << '\n' << __func__ << std::endl;
  auto p = FitLinear<Vect>(xx, uu);
  std::cout << "xx[uu]: ";
  for (size_t i = 0; i < xx.size(); ++i) {
    std::cout << xx[i] << "[" << uu[i] << "]   ";
  }
  std::cout << std::endl;
  std::cout << "g=" << p.first << " u0=" << p.second << std::endl;
}

int main() {
  TestSolveLinear();
  TestFitLinear(
      {Vect(0., 0., 0.), Vect(1., 0., 0.), Vect(0., 1., 0.), Vect(0., 0., 1.)},
      {1., 1., 1., 1.});
  TestFitLinear(
      {Vect(0., 0., 0.), Vect(1., 0., 0.), Vect(0., 1., 0.), Vect(0., 0., 1.)},
      {-1., 1., 1., 1.});
  TestFitLinear(
      {Vect(0., 0., 0.), Vect(1., 0., 0.), Vect(0., 2., 0.), Vect(0., 0., 3.)},
      {0., 1., 1., 1.});
}
