// Created by Petr Karnakov on 21.11.2019
// Copyright 2019 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <iostream>
#include <sstream>
#include <typeinfo>

#include "solver/cond.h"
#include "util/hydro.ipp"

const int dim = 3;
using Scal = double;
using Vect = generic::Vect<Scal, dim>;

std::string P(const void*) {
  return std::string();
}

void TestParse() {
  std::cout << "\nTest Parse" << std::endl;
  BCondAdvection<Scal> cfa;
  cfa.nci = 1;
  ParseAdvectionFaceCond("clear0 1", cfa);
  ParseAdvectionFaceCond("clear1 2", cfa);
  ParseAdvectionFaceCond("halo fill", cfa);
  ParseAdvectionFaceCond("fill_vf 0.4", cfa);
  ParseAdvectionFaceCond("fill_cl 0.6", cfa);
  std::cout << cfa << std::endl;
  ParseAdvectionFaceCond("halo reflect", cfa);
  std::cout << cfa << std::endl;
  for (auto s :
       SplitByDelimiter("clear0 2, clear1 3, fill_vf 0.3, fill_cl 0.7", ',')) {
    ParseAdvectionFaceCond(s, cfa);
  }
  std::cout << cfa << std::endl;
}

int main() {
  TestParse();
}
