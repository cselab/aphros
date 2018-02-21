#include <sstream>
#include <iostream>
#include <cassert>

#include "mesh.hpp"


int main() {
  const int dim = 3;

  using MIdx = geom::MIdxGeneral<dim>;
  using IdxFace = geom::IdxFace;

  const size_t hl = 1;
  MIdx oi(0);
  MIdx si(2);
  MIdx oa= oi - MIdx(hl);
  MIdx sa = si + MIdx(2 * hl);

  geom::BlockFaces<dim> bi(oi, si);
  geom::BlockFaces<dim> ba(oa, sa);

  geom::Range<IdxFace> ra(ba);
  geom::RangeInner<IdxFace, dim> ri(ba, bi);

  for (auto i : ri) {
    std::cout 
      << ba.GetMIdx(i) << " " 
      << ba.GetDirection(i).GetLetter() << " "
      << std::endl;
  }
}
