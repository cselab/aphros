#undef NDEBUG
#include <sstream>
#include <iostream>
#include <cassert>

#include "solver/cond.h"

const int dim = 3;
using Scal = double;
using Vect = GVect<Scal, dim>;

void Test() {
  Rect<Vect> dom{Vect(0), Vect(1)};
  using M = MeshStructured<Scal, dim>;
  MIdx b(0); // lower index
  MIdx s(2);    // size in cells
  int hl = 2;         // halos
  M m = InitUniformMesh<M>(dom, b, s, hl, true, true, s, 0);

  auto bc = m.GetIndexCells();
  auto bf = m.GetIndexFaces();

  auto t = [&](MIdx fw, size_t fd, size_t nci) {
    IdxFace f = bf.GetIdx(fw, Dir(fd));
    IdxCell cmm, cm, cp, cpp;
    GetCellColumn(m, f, nci, cmm, cm, cp, cpp);
    std::cout
        << "fw=" << fw
        << " fd=" << fd
        << " nci=" << nci
        << " cmm=" << bc.GetMIdx(cmm) << " cm=" << bc.GetMIdx(cm)
        << " cp=" << bc.GetMIdx(cp) << " cpp=" << bc.GetMIdx(cpp)
        << std::endl;
  };
  t(MIdx(0), 0, 1);
  t(MIdx(0), 0, 0);
  t(MIdx(0), 1, 1);
  t(MIdx(0), 1, 0);
  t(MIdx(0), 2, 1);
  t(MIdx(0), 2, 0);
  t(MIdx(1), 0, 1);
  t(MIdx(1), 0, 0);
}

int main() {
  Test();
}
