// Created by Petr Karnakov on 07.02.2020
// Copyright 2020 ETH Zurich

#undef NDEBUG
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <random>
#include <sstream>

#include "geom/vect.h"
#include "solver/approx_eb.h"
#include "solver/embed.h"

using Scal = double;
using Vect = generic::Vect<Scal, 3>;

void TestSolveLinear() {
  constexpr size_t N = 4;
  std::array<Scal, N* N> a = {
      0, 0, 0.1, 9.8, 0.7, 1, 0.8, 0.3, 1, 0, 0.3, 0.2, 1, 0, 0.3, 0.1,
  };
  std::array<Scal, N> b = {
      1,
      2,
      8,
  };
  auto x = ULinearFit<Vect>::SolveLinear(a, b);
  using V = generic::Vect<Scal, N>;
  std::cout << "x=" << V(x) << std::endl;
  std::cout << "b=" << V(b) << std::endl;
  std::cout << "a*x=" << V(ULinearFit<Vect>::Mul(a, x)) << std::endl;
}

template <class T>
void TestFitLinear(const std::vector<Vect>& xx, const std::vector<T>& uu) {
  std::cout << '\n' << __func__ << std::endl;
  auto p = ULinearFit<Vect>::FitLinear(xx, uu);
  std::cout << "xx[uu]: ";
  for (size_t i = 0; i < xx.size(); ++i) {
    std::cout << xx[i] << "[" << uu[i] << "]   ";
  }
  std::cout << std::endl;
  std::cout << "g=" << p.first << " u0=" << p.second << std::endl;
}

int main() {
  TestSolveLinear();
  TestFitLinear<Scal>(
      {Vect(0., 0., 0.), Vect(1., 0., 0.), Vect(0., 1., 0.), Vect(0., 0., 1.)},
      {1., 1., 1., 1.});
  TestFitLinear<Scal>(
      {Vect(0., 0., 0.), Vect(1., 0., 0.), Vect(0., 1., 0.), Vect(0., 0., 1.)},
      {-1., 1., 1., 1.});
  TestFitLinear<Scal>(
      {Vect(0., 0., 0.), Vect(1., 0., 0.), Vect(0., 2., 0.), Vect(0., 0., 3.)},
      {0., 1., 1., 1.});

  TestFitLinear<Vect>(
      {Vect(0., 0., 0.), Vect(1., 0., 0.), Vect(0., 2., 0.), Vect(0., 0., 3.)},
      {Vect(0., 0., 0.), Vect(1., 2., 3.), Vect(1., 2., 3.), Vect(1., 2., 3.)});
}
