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
  Vect m(1., 1., 1.);  // normal 
  m /= m.norm();

  Scal c;   // volume fraction
  Scal area; // area 
  Scal alpha; // line constant (length measured from corner)

  c = 0.5;
  alpha = gfs_plane_alpha(m, c);
  area = gfs_plane_volume(m, alpha);
  Pr(m, c, area, alpha);

  return 0;
}
