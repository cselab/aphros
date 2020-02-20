// Created by Petr Karnakov on 20.02.2020
// Copyright 2020 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>

#include "parse/codeblocks.h"

void Test() {
  std::cout << "\n" << __func__ << std::endl;

  auto bb = ParseCodeBlocks(std::cin);

  size_t i = 0;
  for (auto b : bb) {
    std::cout << "i=" << i << std::endl;
    std::cout << "name='" << b.name << "'" << std::endl;
    std::cout << "content='" << b.content << "'" << std::endl;
    std::cout << std::endl;
    ++i;
  }
}

int main() {
  Test();
}
