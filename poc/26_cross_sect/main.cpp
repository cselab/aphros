// Created by Petr Karnakov on 07.01.2020
// Copyright 2020 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <iostream>
#include <memory>
#include <string>

#include <geom/vect.h>
#include <solver/reconst.h>

using Scal = double;
using Vect = GVect<Scal, 3>;
using R = Reconst<Scal>;

int main() {
  Vect x0, x1;
  bool q = R::PolyInter(
      Vect{0., 0., 0.}, Vect{1., 1., 1.}, 0.3, Vect{1., 1., 1.},
      Vect{0., 0., 0.}, Vect{0., 1., 0.}, x0, x1);
  std::cout << x0 << std::endl;
  std::cout << x1 << std::endl;
  std::cout << q << std::endl;
}
