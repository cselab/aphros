#include <cassert>
#include <iostream>
#include <sstream>
#include <fstream>

#include "parse/vars.h"
#include "parse/interp.h"

namespace simple {

void Simple() {
  Vars par;
  Interp ip(par);

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
  Interp ip(par);
  ip.RunAll(fi);

  // print variables to stringstream
  std::stringstream so;
  ip.PrintAll(so);

  // echo to stderr
  std::cerr << so.str();

  // read reference data
  std::string nr = "a.out";
  std::ifstream fr(nr);
  std::stringstream sr;
  sr << fr.rdbuf();

  // compare
  assert(so.str() == sr.str());
}

int main() {
  simple::Simple();

  TestFile();
}
