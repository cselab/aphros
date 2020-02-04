// Created by Petr Karnakov on 04.02.2020
// Copyright 2020 ETH Zurich

#include <iostream>

#include <util/sysinfo.h>

int main() {
  std::cout << sysinfo::GetMem() / double((1 << 20)) << std::endl;
}
