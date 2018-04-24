#pragma once

#include <string>
#include <fstream>

#include "hydro/block.h"
#include "hydro/field.h"

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

  auto s = b.GetDimensions();
  o << s[0] << " " << s[1] << " " << s[2] << std::endl;

  for (auto c : b.Range()) {
    o << u[c] << " ";
  }

  o << std::endl;
}
