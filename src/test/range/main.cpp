// Created by Petr Karnakov on 01.01.2020
// Copyright 2020 ETH Zurich

#undef NDEBUG
#include <stdlib.h>
#include <cassert>
#include <iostream>
#include <map>

#include "geom/filter.h"
#include "geom/range.h"
#include "geom/transform.h"

#define E(x)                   \
  {                            \
    std::cout << (#x) << "\n"; \
    x;                         \
    std::cout << "\n\n";       \
  }

int main() {
  E(for (auto i : GRange<size_t>(10, 20)) { std::cout << i << " "; });
  E(for (auto i
         : MakeFilterIterator(GRange<size_t>(10, 20), [](size_t a) {
           return a % 2 == 0;
         })) { std::cout << i << " "; });
  E(for (auto i
         : MakeFilterIterator(
             MakeFilterIterator(
                 GRange<size_t>(10, 20), [](size_t a) { return a % 2 == 0; }),
             [](size_t a) { return a % 3 == 0; })) { std::cout << i << " "; });
  E(for (auto i
         : MakeTransformIterator<size_t>(GRange<size_t>(10, 20), [](size_t a) {
           return a;
         })) { std::cout << i << " "; });
  E(for (auto i
         : MakeTransformIterator<double>(GRange<size_t>(10, 20), [](size_t a) {
           return a + 0.5;
         })) { std::cout << i << " "; });
}
