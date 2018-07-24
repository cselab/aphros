#pragma once

#include <string>
#include <fstream>
#include <vector>

#include "geom/block.h"
#include "geom/field.h"
#include "geom/rangein.h"

// Dump values to text file. 
// u: scalar field defined on b
// b: block of cells
// op: output path
// Format:
// <Nx> <Ny> <Nz>
// <data:x=0,y=0,z=0> <data:x=1,y=0,z=0> ...
template <class Scal>
void Dump(const FieldCell<Scal>& u, 
          const GBlock<IdxCell, 3>& b,
          std::string op) {
  std::ofstream o;
  o.open(op);
  o.precision(20);

  auto s = b.GetSize();
  o << s[0] << " " << s[1] << " " << s[2] << std::endl;

  for (auto c : b.Range()) {
    o << u[c] << " ";
  }

  o << std::endl;
}

