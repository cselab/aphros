#include <sstream>
#include <iostream>
#include <cassert>

#include "mesh.hpp"
#include "mesh3d.hpp"

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

const int dim = 3;
using MIdx = geom::GMIdx<dim>;
using IdxFace = geom::IdxFace;
using Dir = geom::GDir<dim>;
using Scal = double;
using Vect = geom::GVect<Scal, dim>;

void TestBlock() {
  const size_t hl = 1;
  MIdx oi(0); // origin inner
  MIdx si(2); // size inner
  MIdx oa= oi - MIdx(hl); // origin all
  MIdx sa = si + MIdx(2 * hl); // size all

  geom::GBlockFaces<dim> bi(oi, si);
  geom::GBlockFaces<dim> ba(oa, sa);

  geom::GRange<IdxFace> ra(ba);
  geom::GRangeIn<IdxFace, dim> ri(ba, bi);

  const MIdx xp0 = oa - Dir(0);
  MIdx xp = xp0;
  Dir dp(0); // direction

  // Check that whole inner block covered with ascending indices
  for (auto i : ri) {
    auto x = ba.GetMIdx(i);
    auto d = ba.GetDir(i);

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

void TestMesh() {
  geom::Rect<Vect> dom(Vect(0), Vect(1));
  using M = geom::MeshStructured<Scal, dim>;
  M m = geom::InitUniformMesh<M>(dom, MIdx(0), MIdx(5), 1);

  //for (auto i : m.GetAll<geom::IdxFace>()) {
  int a = 0;
  for (auto i : m.Faces()) {
    std::cout << i.GetRaw() << std::endl;
    ++a;
  }
  std::cout << a << std::endl;

}

int main() {
  TestBlock();

  TestMesh();
}
