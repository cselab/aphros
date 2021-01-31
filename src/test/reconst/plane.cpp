// Created by Petr Karnakov on 01.09.2019
// Copyright 2019 ETH Zurich

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>

#include "geom/vect.h"
#include "solver/reconst.h"

using Scal = double;
using Vect = generic::Vect<Scal, 3>;
using Vect2 = generic::Vect<Scal, 2>;
using R = Reconst<Scal>;

int main(int, char** argc) {
  std::string f = argc[1];
  Scal nx, ny, nz;
  Scal v;
  nx = atof(argc[2]);
  ny = atof(argc[3]);
  nz = atof(argc[4]);
  v = atof(argc[5]);

  Scal r = 0;
  if (f == "alpha_ch") {
    r = R::GetLineA1(Vect(nx, ny, nz), v);
  } else if (f == "volume_ch") {
    r = R::GetLineU1(Vect(nx, ny, nz), v);
  }
  std::cout << r << std::endl;
}
