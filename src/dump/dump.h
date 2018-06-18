#pragma once

#include <string>
#include <fstream>

#include "geom/block.h"
#include "geom/field.h"

// Dump values in inner cells to text file. 
// u: scalar field
// b: block
// op: output path
// Format:
// <Nx> <Ny> <Nz>
// <data:x=0,y=0,z=0> <data:x=1,y=0,z=0> ...
template <class Scal, class B=GBlock<IdxCell, 3>>
void Dump(const FieldCell<Scal>& u, const B& b, std::string op) {
  std::ofstream o(op.c_str());

  auto l = o.flags();
  o.precision(20);

  auto s = b.GetDimensions();
  o << s[0] << " " << s[1] << " " << s[2] << std::endl;

  for (auto c : b.Range()) {
    o << u[c] << " ";
  }

  o.flags(l);

  o << std::endl;
}
