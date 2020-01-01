#undef NDEBUG
#include <stdlib.h>
#include <cassert>
#include <iostream>
#include <map>

#include "geom/range.h"
#include "geom/filter.h"

#define E(x) { std::cout << (#x) << std::endl; (x); std::cout << std::endl; }

int main(int an, char* av[]) {
  E(for (auto i : GRange<size_t>(10, 20)) {
      std::cout << i << " ";
  });
}
