// Created by Petr Karnakov on 17.07.2018
// Copyright 2018 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <iostream>

#include "util/logger.h"
#include "dump/dumper.h"

void TestName() {
  std::cout << NAMEVALUE(GetDumpName("vf", ".h5", 8, 100)) << '\n';
  std::cout << NAMEVALUE(GetDumpName("vf", ".h5", 8)) << '\n';
}

int main() {
  TestName();
}
