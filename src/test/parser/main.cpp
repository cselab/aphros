// Created by Petr Karnakov on 30.05.2018
// Copyright 2018 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>

#include "parse/parser.h"
#include "parse/vars.h"

namespace simple {

void Simple() {
  std::cout << "\n" << __func__ << std::endl;
  Vars par;
  Parser ip(par);

  std::stringstream s;
  s << "set string a 1" << std::endl;
  s << "set int b 1" << std::endl;
  s << "set double c 1" << std::endl;
  s << "set vect d 1" << std::endl;

  ip.RunAll(s);
  ip.PrintAll(std::cout);

  assert(par.String["a"] == "1");
  assert(par.Int["b"] == 1);
  assert(par.Double["c"] == 1.);
  assert(par.Vect["d"] == std::vector<double>({1}));
}

} // namespace simple

void TestFile() {
  std::cout << "\n" << __func__ << std::endl;
  // open script
  std::string ni = "a.conf";
  std::ifstream fi(ni);

  // run all commands
  Vars par;
  Parser ip(par);
  ip.RunAll(fi);

  // print variables to cout
  ip.PrintAll(std::cout);
}

int main() {
  simple::Simple();

  TestFile();
}
