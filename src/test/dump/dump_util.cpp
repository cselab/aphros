// Created by Petr Karnakov on 17.07.2018
// Copyright 2018 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <iostream>
#include <fstream>

#include "dump/dumper.h"
#include "dump/raw.h"
#include "dump/raw.ipp"
#include "util/logger.h"
#include "geom/vect.h"

void TestName() {
  std::cout << NAMEVALUE(GetDumpName("vf", ".h5", 8, 100)) << '\n';
  std::cout << NAMEVALUE(GetDumpName("vf", ".h5", 8)) << '\n';
}

struct M {
  static constexpr size_t dim = 3;
  using Scal = double;
  using Vect = generic::Vect<Scal, dim>;
  using MIdx = generic::Vect<int, dim>;
};

void TestReadXmf() {
  const std::string path = "input.xmf";
  std::ifstream fin(path);
  using Raw = dump::Raw<M>;
  const auto meta = Raw::ReadXmf(fin);

  std::cout << NAMEVALUE(meta.name) << std::endl;
  std::cout << NAMEVALUE(meta.binpath) << std::endl;
  std::cout << NAMEVALUE(Raw::TypeToString(meta.type)) << std::endl;
  std::cout << NAMEVALUE(meta.dimensions) << std::endl;
  std::cout << NAMEVALUE(meta.start) << std::endl;
  std::cout << NAMEVALUE(meta.stride) << std::endl;
  std::cout << NAMEVALUE(meta.count) << std::endl;
  std::cout << NAMEVALUE(meta.seek) << std::endl;
  std::cout << NAMEVALUE(meta.origin) << std::endl;
  std::cout << NAMEVALUE(meta.spacing) << std::endl;

  dump::Raw<M>::WriteXmf(std::cout, meta);
}

int main() {
  TestName();
  TestReadXmf();
}
