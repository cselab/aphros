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

int main() {
  using Scal = double;
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
