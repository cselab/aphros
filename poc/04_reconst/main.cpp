#include <iostream>
#include <cstdlib>
#include <cmath>

#include "reconst.h"

void Pr(Vect m, Scal c, Scal area, Scal alpha) {
  std::cout
      << " m=" << m
      << " c=" << c
      << " area=" << area
      << " alpha=" << alpha
      << std::endl;
}

int main() {
  Vect m(1., 0., 0.);  // normal 
  m /= m.norm();

  Scal c;   // volume fraction
  Scal area; // area 
  Scal alpha; // line constant (length measured from corner)

  c = 0.5;
  alpha = gfs_line_alpha(m, c);
  area = gfs_line_area(m, alpha);
  Pr(m, c, area, alpha);

  alpha = std::sqrt(2) * 0.5;
  c = gfs_line_area(m, alpha);
  area = 0;
  Pr(m, c, area, alpha);

  return 0;
}
