#undef NDEBUG
#include <stdlib.h>
#include <cassert>
#include <iostream>
#include <map>

#include "geom/filter.h"
#include "geom/range.h"

#define E(x)                   \
  {                            \
    std::cout << (#x) << "\n"; \
    x;                         \
    std::cout << "\n\n";       \
  }

int main() {
  E(for (auto i : GRange<size_t>(10, 20)) { std::cout << i << " "; });
  E(for (auto i
         : MakeFilter(GRange<size_t>(10, 20), [](size_t a) {
           return a % 2 == 0;
         })) { std::cout << i << " "; });
  E(for (auto i : MakeFilter(MakeFilter(GRange<size_t>(10, 20), [](size_t a){
    return a % 2 == 0; }), [](size_t a){
    return a % 3 == 0; })) {
    std::cout << i << " ";
  });
}
