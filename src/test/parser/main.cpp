// Created by Petr Karnakov on 30.05.2018
// Copyright 2018 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>

#include "parse/parser.h"
#include "parse/vars.h"
#include "parse/config.h"

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

struct Config : public ConfigBase {
  VAR_DOUBLE(height);
  VAR_INT(size);
  VAR_VECT(elems);
  VAR_VECT3(gravity);
  VAR_STRING(name);
  VAR_BOOL(enable_fluid);
};


void TestConfig() {
  std::cout << "\n" << __func__ << std::endl;
  Vars var;
  Parser parser(var);

  std::stringstream s;
  s << R"EOF(

set double height 1.2
set int size 3
set vect elems 1 2 3
set vect gravity 4 5 6
set string name name
set int enable_fluid 3

)EOF" << std::endl;

  parser.RunAll(s);

  using Vect = generic::Vect<double, 3>;

  Config config;
  config.Read(var);
  std::cout << config.height << std::endl;
  std::cout << config.size << std::endl;
  std::cout << Vect(config.elems) << std::endl;
  std::cout << config.gravity << std::endl;
  std::cout << config.name << std::endl;
  std::cout << config.enable_fluid << std::endl;
}

int main() {
  simple::Simple();

  TestFile();
  TestConfig();
}
