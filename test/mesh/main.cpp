#include <sstream>
#include <iostream>
#include <cassert>

#include "mesh.hpp"

// Returns true if a < b (lex starting from end)
template <class Vect>
bool Cmp(Vect a, Vect b) {
  int i = Vect::dim;
  while (i--) {
    if (a[i] != b[i]) {
      return a[i] < b[i];
    }
  }
  return false;
}

int main() {
  const int dim = 3;

  using MIdx = geom::MIdxGeneral<dim>;
  using IdxFace = geom::IdxFace;
  using Dir = geom::Direction<dim>;

  const size_t hl = 1;
  MIdx oi(0); // origin inner
  MIdx si(2); // size inner
  MIdx oa= oi - MIdx(hl); // origin all
  MIdx sa = si + MIdx(2 * hl); // size all

  geom::BlockFaces<dim> bi(oi, si);
  geom::BlockFaces<dim> ba(oa, sa);

  geom::Range<IdxFace> ra(ba);
  geom::RangeInner<IdxFace, dim> ri(ba, bi);

  const MIdx xp0 = oa - Dir(0);
  MIdx xp = xp0;
  Dir dp(0); // direction

  // Check that whole inner block covered with ascending indices
  for (auto i : ri) {
    auto x = ba.GetMIdx(i);
    auto d = ba.GetDirection(i);

    // Next direction, reset xp
    if (dp < d) {
      xp = xp0;
    }

    // std::cerr << x << " " << d.GetLetter() << std::endl;
    assert(Cmp(xp, x));
    assert(oi <= x && x < oi + si + d);

    xp = x;
    dp = d;
  }
}
