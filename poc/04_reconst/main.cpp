#include <iostream>
#include <cstdlib>

#include "reconst.h"

void Pr(Vect m, Scal a, Scal s) {
  std::cout << "m=" << m << " a=" << a << " s=" << s << std::endl;
}

int main(int c, char** v) {
  Vect m(1., 1., 0.);  // normal 
  m /= m.norm();
  Scal a = 0.5;   // volume fraction
  Scal s; // area

  s = gfs_line_area(m, a);
  Pr(m, a, s);

  return 0;
}
