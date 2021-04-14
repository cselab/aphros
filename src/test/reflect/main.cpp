// Created by Petr Karnakov on 11.11.2019
// Copyright 2019 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <iostream>
#include <sstream>

#include "geom/mesh.h"
#include "solver/solver.h"

const int dim = 3;
using MIdx = generic::MIdx<dim>;
using Dir = GDir<dim>;
using Scal = double;
using Vect = generic::Vect<Scal, dim>;

#define CMP(a, b) assert(Cmp(a, b));

// Print CMP
#define PCMP(a, b)                                                    \
  std::cerr << #a << "=" << a << ", " << #b << "=" << b << std::endl; \
  CMP(a, b);

void TestCellColumn() {
  Rect<Vect> dom{Vect(0), Vect(1)};
  using M = MeshCartesian<Scal, dim>;
  MIdx b(0); // lower index
  MIdx s(2); // size in cells
  int hl = 2; // halos
  const M m{b, s, dom, hl, true, true, s, 0};

  auto bc = m.GetIndexCells();
  auto bf = m.GetIndexFaces();

  auto t = [&](MIdx fw, size_t fd, size_t nci) {
    IdxFace f = bf.GetIdx(fw, Dir(fd));
    IdxCell cmm, cm, cp, cpp;
    GetCellColumn(m, f, nci, cmm, cm, cp, cpp);
    std::cout << "fw=" << fw << " fd=" << fd << " nci=" << nci
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
  TestCellColumn();
}
