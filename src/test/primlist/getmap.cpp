// Created by Petr Karnakov on 22.09.2019
// Copyright 2019 ETH Zurich

#undef NDEBUG
#include <stdlib.h>
#include <cassert>
#include <iostream>
#include <map>

#include "func/primlist.ipp"

using Scal = double;
using Vect = generic::Vect<Scal, 3>;

int main(int an, char* av[]) {
  if (an < 5) {
    std::cerr << "usage: " << av[0] << " name values keys nrequired"
              << std::endl;
    return 1;
  }
  try {
    std::map<std::string, Scal> r =
        Imp<Vect>::GetMap(av[1], av[2], av[3], atoi(av[4]));
    for (auto k : r) {
      std::cout << k.first << ":" << k.second << " ";
    }
    std::cout << std::endl;
  } catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
