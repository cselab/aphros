#undef NDEBUG
#include <cassert>
#include <iostream>
#include <map>
#include <stdlib.h>

#include "func/primlist.h"

using Scal = double;

int main(int an, char* av[]) {
  std::map<std::string, Scal> r =
      UPrimList<Scal>::Parse(av[1], av[2], av[3], atoi(av[4]));
  for (auto k : r) {
    std::cout << k.first << ":" << k.second << " ";
  }
  std::cout << std::endl;
}
