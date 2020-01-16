// Created by Petr Karnakov on 20.07.2018
// Copyright 2018 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

#include "util/sysinfo.h"

using namespace sysinfo;

volatile char g;

// Perform action allocating s bytes of memory
void F(size_t s) {
  static std::vector<char> v(s, 0);
}

double MB(size_t s) {
  return s / (double(1 << 20));
}

void Test() {
  size_t a, b;
  size_t s = (1 << 20) * 100; // memory to allocate

  a = GetMem();

  F(s);

  b = GetMem();
  std::cout << MB(a) << " " << MB(b) << std::endl;
  assert(b >= a + s);
}

int main() {
  Test();
}
