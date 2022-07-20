// Created by Petr Karnakov on 17.07.2018
// Copyright 2018 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <fstream>
#include <iostream>

#include "dump/dumper.h"
#include "dump/xmf.h"
#include "dump/xmf.ipp"
#include "geom/vect.h"
#include "util/logger.h"

void TestName() {
  std::cout << NAMEVALUE(GetDumpName("vf", ".h5", 8, 100)) << '\n';
  std::cout << NAMEVALUE(GetDumpName("vf", ".h5", 8)) << '\n';
}

static constexpr size_t dim = 3;
using Scal = double;
using Vect = generic::Vect<Scal, dim>;

void TestReadXmf() {
  const std::string path = "input.xmf";
  std::ifstream fin(path);
  using Xmf = dump::Xmf<Vect>;
  const auto meta = Xmf::ReadXmf(fin);

  std::cout << NAMEVALUE(meta.name) << std::endl;
  std::cout << NAMEVALUE(meta.binpath) << std::endl;
  std::cout << NAMEVALUE(dump::TypeToString(meta.type)) << std::endl;
  std::cout << NAMEVALUE(meta.dimensions) << std::endl;
  std::cout << NAMEVALUE(meta.start) << std::endl;
  std::cout << NAMEVALUE(meta.stride) << std::endl;
  std::cout << NAMEVALUE(meta.count) << std::endl;
  std::cout << NAMEVALUE(meta.seek) << std::endl;
  std::cout << NAMEVALUE(meta.origin) << std::endl;
  std::cout << NAMEVALUE(meta.spacing) << std::endl;

  Xmf::WriteXmf(std::cout, meta, false);
}

int main() {
  TestName();
  TestReadXmf();
}
