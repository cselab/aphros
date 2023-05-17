// Created by Petr Karnakov on 01.09.2019
// Copyright 2019 ETH Zurich

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>

#include "func/init_u.h"
#include "func/init.ipp"
#include "geom/vect.h"

using Scal = double;
using Vect = generic::Vect<Scal, 3>;

int main(int, char** argc) {
  std::string f = argc[1];
  if (f == "area") {
    const std::array<Scal, 4> ee{
        atof(argc[2]), atof(argc[3]), atof(argc[4]), atof(argc[5])};
    std::cout << GetFaceAreaFraction(ee);
  }
  std::cout << std::endl;
}
