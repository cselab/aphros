// Created by Petr Karnakov on 03.02.2020
// Copyright 2020 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <random>
#include <sstream>

#include "dump/dump.h"
#include "geom/mesh.h"
#include "solver/approx.h"
#include "solver/approx_eb.h"
#include "solver/cond.h"
#include "solver/embed.h"
#include "solver/solver.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;
auto Dirs = GRange<size_t>(3);
using EB = Embed<M>;
using UEB = UEmbed<M>;

// h: cell size
std::unique_ptr<M> CreateMesh() {
  const MIdx size(16);
  Rect<Vect> dom(Vect(0), Vect(1));
  return std::make_unique<M>(
      InitUniformMesh<M>(dom, MIdx(0), size, 2, true, true, size, 0));
}

std::unique_ptr<EB> CreateEmbed(M& m) {
  auto peb = std::make_unique<EB>(m, 0);
  FieldNode<Scal> fnl(m);
  auto block = m.GetInBlockCells().GetSize();
  auto h = m.GetCellSize();
  for (auto n : m.AllNodes()) {
    const auto x = m.GetNode(n) / h / Vect(block);
    auto dx = Vect(0.5) - x;
    fnl[n] = 0.41 - dx.norm();
  }
  do {
    peb->Init(fnl);
  } while (m.Pending());
  return peb;
}

void Test() {
  auto pm = CreateMesh();
  auto& m = *pm;
}

int main() {
  Test();
}
