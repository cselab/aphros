// Created by Petr Karnakov on 22.09.2019
// Copyright 2019 ETH Zurich

#undef NDEBUG
#include <stdlib.h>
#include <cassert>
#include <iostream>
#include <map>
#include <sstream>

#include "func/primlist.h"

using Vect = generic::Vect<double, 3>;

int main(int an, char* av[]) {
  if (an < 2) {
    std::cerr << "usage: " << av[0] << " primitive" << std::endl;
    return 1;
  }
  using U = UPrimList<Vect>;
  try {
    std::stringstream ss(av[1]);
    auto pp = U::GetPrimitives(ss, 3);
    if (pp.size() == 1) {
      std::cout << pp[0] << std::endl;
    } else {
      std::cout << "pp.size()=" << pp.size() << std::endl;
    }
  } catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
