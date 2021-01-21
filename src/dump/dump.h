// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <fstream>
#include <string>
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
void Dump(
    const FieldCell<Scal>& u, const GIndex<IdxCell, 3>& ndc,
    const GBlock<IdxCell, 3>& bc, std::string op) {
  std::ofstream o;
  o.open(op);
  o.precision(20);

  auto s = bc.GetSize();
  o << s[0] << " " << s[1] << " " << s[2] << std::endl;

  for (auto c : GRangeIn<IdxCell, 3>(ndc, bc)) {
    o << u[c] << " ";
  }

  o << std::endl;
}

// Creates csv file with scalar fields gathered from all blocks.
// data: pairs {name, values}
// path: output path
template <class Scal>
void DumpCsv(
    const std::vector<std::pair<std::string, std::vector<Scal>>>& data,
    std::string path);

// Creates csv file with scalar fields gathered from all blocks.
// data: pairs {name, values}
// path: output path
template <class M>
void DumpCsv(
    const std::vector<std::pair<std::string, std::vector<typename M::Scal>>>&
        data,
    std::string path, M& m);
