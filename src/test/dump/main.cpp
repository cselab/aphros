// Created by Petr Karnakov on 17.07.2018
// Copyright 2018 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <iostream>

#include "dump/dumper.h"

void TestName() {
  assert(GetDumpName("vf", ".h5", 8, 100) == "vf_0008_0100.h5");
  assert(GetDumpName("vf", ".h5", 8) == "vf_0008.h5");
}

int main() {
  TestName();
}
