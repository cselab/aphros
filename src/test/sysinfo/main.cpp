// Created by Petr Karnakov on 20.07.2018
// Copyright 2018 ETH Zurich

#undef NDEBUG
#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

#include "util/filesystem.h"
#include "util/format.h"
#include "util/logger.h"
#include "util/sysinfo.h"

using namespace sysinfo;

volatile char g;

// Perform action allocating size bytes of memory
void F(size_t size) {
  static std::vector<char> v(size, 0);
}

double MiB(size_t size) {
  return size / (double(1 << 20));
}

void Test() {
  size_t a, b;
  size_t size = (1 << 20) * 100; // memory to allocate

  a = GetMem();

  F(size);

  b = GetMem();
  std::cerr << util::Format(
      "Allocated memory (MiB): {:.0f}\nMeasured increase (MiB): {:.0f}\n",
      MiB(size), MiB(b - a));
  const auto factor = 0.5;
  fassert(
      b - a >= size * factor, util::Format(
                                  "required at least {:.0f}, got {:.0f}",
                                  MiB(factor * size), MiB(b - a)));
}

std::ostream& operator<<(
    std::ostream& out, const std::array<std::string, 2>& v) {
  out << "(" << v[0] << ", " << v[1] << ")";
  return out;
}

void TestFilesystem() {
  using namespace util;
  std::cout << NAMEVALUE(GetBasename("")) << std::endl;
  std::cout << NAMEVALUE(GetBasename("./")) << std::endl;
  std::cout << NAMEVALUE(GetBasename("./.")) << std::endl;
  std::cout << NAMEVALUE(GetBasename("a/b")) << std::endl;
  std::cout << NAMEVALUE(GetBasename("a")) << std::endl;
  std::cout << NAMEVALUE(GetDirname("")) << std::endl;
  std::cout << NAMEVALUE(GetDirname("./")) << std::endl;
  std::cout << NAMEVALUE(GetDirname("./.")) << std::endl;
  std::cout << NAMEVALUE(GetDirname("a/b")) << std::endl;
  std::cout << NAMEVALUE(GetDirname("a")) << std::endl;
  std::cout << NAMEVALUE(IsFile("test")) << std::endl;
  std::cout << NAMEVALUE(IsFile("main2.cpp")) << std::endl;
  std::cerr << NAMEVALUE(GetRealpath("test")) << std::endl;
  std::cout << NAMEVALUE(GetRealpath("main2.cpp")) << std::endl;
  std::cout << NAMEVALUE(Join("a", "")) << std::endl;
  std::cout << NAMEVALUE(Join("", "b")) << std::endl;
  std::cout << NAMEVALUE(Join("a", "b")) << std::endl;
  std::cout << NAMEVALUE(Join("a/", "b")) << std::endl;
  std::cout << NAMEVALUE(Join("a", "/b")) << std::endl;
  std::cout << NAMEVALUE(IsDir(".")) << std::endl;
  std::cout << NAMEVALUE(SplitExt(".a")) << std::endl;
  std::cout << NAMEVALUE(SplitExt("a.b")) << std::endl;
  std::cout << NAMEVALUE(SplitExt("a/.b")) << std::endl;
  std::cout << NAMEVALUE(SplitExt("a/b.c")) << std::endl;
  std::cout << NAMEVALUE(SplitExt("a/b.c.d")) << std::endl;
  std::cout << NAMEVALUE(SplitExt("a")) << std::endl;
}

int main() {
  TestFilesystem();
  Test();
}
