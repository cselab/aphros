// Created by Petr Karnakov on 20.07.2018
// Copyright 2018 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

#include "util/sysinfo.h"
#include "util/filesystem.h"
#include "util/logger.h"

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
  std::cerr << MB(a) << " " << MB(b) << std::endl;
  assert(b >= a + s);
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
  std::cout << NAMEVALUE(IsFile("main.cpp")) << std::endl;
  std::cout << NAMEVALUE(IsFile("main2.cpp")) << std::endl;
  std::cerr << NAMEVALUE(GetRealpath("main.cpp")) << std::endl;
  std::cout << NAMEVALUE(GetRealpath("main2.cpp")) << std::endl;
  std::cout << NAMEVALUE(Join("a", "")) << std::endl;
  std::cout << NAMEVALUE(Join("", "b")) << std::endl;
  std::cout << NAMEVALUE(Join("a", "b")) << std::endl;
  std::cout << NAMEVALUE(Join("a/", "b")) << std::endl;
  std::cout << NAMEVALUE(Join("a", "/b")) << std::endl;
  std::cout << NAMEVALUE(IsDir(".")) << std::endl;
}

int main() {
  TestFilesystem();
  Test();
}
