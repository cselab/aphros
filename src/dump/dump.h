#pragma once

#include <string>
#include <fstream>
#include <vector>

#include "geom/block.h"
#include "geom/field.h"
#include "geom/rangein.h"

// Dump values to text file. 
// u: scalar field defined on b
// ndc: index cells
// bc: block cells
// op: output path
// Format:
// <Nx> <Ny> <Nz>
// <data:x=0,y=0,z=0> <data:x=1,y=0,z=0> ...
template <class Scal>
void Dump(const FieldCell<Scal>& u, const GIndex<IdxCell, 3>& ndc,
          const GBlock<IdxCell, 3>& bc, std::string op) {
  std::ofstream o;
  o.open(op);
  o.precision(20);
  o.sync_with_stdio(false);

  auto s = bc.GetSize();
  o << s[0] << " " << s[1] << " " << s[2] << std::endl;

  for (auto c : GRangeIn<IdxCell, 3>(ndc, bc) ) {
    o << u[c] << " ";
  }

  o << std::endl;
}

