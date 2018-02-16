#include <cassert>
#include <iostream>
#include <sstream>

#include "Vars.h"
#include "Interp.h"

namespace simple {

void Simple() {
  Vars par;
  Interp interp(par);

  std::stringstream s;
  s << "set string a 1" << std::endl;
  s << "set int b 1" << std::endl;
  s << "set double c 1" << std::endl;
  s << "set vect d 1" << std::endl;

  interp.RunAll(s);

  assert(par.String["a"] == "1");
  assert(par.Int["b"] == 1);
  assert(par.Double["c"] == 1.);
  assert(par.Vect["d"] == std::vector<double>({1}));
}

} // namespace simple


int main() {
  simple::Simple();
}
