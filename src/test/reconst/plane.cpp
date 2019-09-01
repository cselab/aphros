#include <iostream>
#include <sstream>
#include <string>
#include <stdio.h>

#include "geom/vect.h"
#include "solver/reconst.h"

using Scal = double;
using Vect = GVect<Scal, 3>;
using Vect2 = GVect<Scal, 2>;
using R = Reconst<Scal, false>;

int main(int, char** argc) {
  std::string f = argc[1];
  Scal nx, ny, nz;
  Scal v;
  nx = atof(argc[2]);
  ny = atof(argc[3]);
  nz = atof(argc[4]);
  v = atof(argc[5]);

  Scal r;
  if (f == "alpha_ba") {
    r = R::plane_alpha(v, nx, ny, nz);
  } else if (f == "alpha_ch") {
    r = R::GetLineA1(Vect(nx, ny, nz), v);
  } else if (f == "volume_ba") {
    r = R::plane_volume(nx, ny, nz, v);
  } else if (f == "volume_ch") {
    r = R::GetLineU1(Vect(nx, ny, nz), v);
  }
  std::cout << r << std::endl;
}
