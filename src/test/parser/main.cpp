#undef NDEBUG
#include <cassert>
#include <iostream>
#include <sstream>
#include <fstream>

#include "parse/vars.h"
#include "parse/parser.h"

namespace simple {

void Simple() {
  Vars par;
  Parser ip(par);

  std::stringstream s;
  s << "set string a 1" << std::endl;
  s << "set int b 1" << std::endl;
  s << "set double c 1" << std::endl;
  s << "set vect d 1" << std::endl;

  ip.RunAll(s);
  ip.PrintAll();

  assert(par.String["a"] == "1");
  assert(par.Int["b"] == 1);
  assert(par.Double["c"] == 1.);
  assert(par.Vect["d"] == std::vector<double>({1}));
}

} // namespace simple


void TestFile() {
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
