// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <fstream>
#include <string>
#include <vector>

#include "geom/block.h"
#include "geom/field.h"
#include "geom/rangein.h"
#include "parse/util.h"

namespace dump {

// Dump values to text file.
// u: scalar field defined on b
// ndc: index cells
// bc: block cells
// op: output path
// Format:
// <Nx> <Ny> <Nz>
// <data:x=0,y=0,z=0> <data:x=1,y=0,z=0> ...
template <class Scal, size_t dim>
void DumpFieldPlain(
    const FieldCell<Scal>& u, const GIndex<IdxCell, dim>& ndc,
    const GBlock<IdxCell, dim>& bc, std::string op) {
  std::ofstream o;
  o.open(op);
  o.precision(20);

  o << bc.GetSize().to_string() << std::endl;

  for (auto c : GRangeIn<IdxCell, dim>(ndc, bc)) {
    o << u[c] << " ";
  }

  o << std::endl;
}

// Reads CSV file. Returns columns as pairs {name, value}.
// path: output path
template <class Scal>
std::vector<std::pair<std::string, std::vector<Scal>>> ReadCsv(
    std::istream& in, char delim = ',');
template <class Scal>
std::vector<std::pair<std::string, std::vector<Scal>>> ReadCsv(
    const std::string& path, char delim = ',');

// Creates CSV file.
// data: pairs {name, values}
// path: output path
template <class Scal>
void DumpCsv(
    const std::vector<std::pair<std::string, std::vector<Scal>>>& data,
    std::string path, char delim=',');

// Creates one CSV file with scalar fields gathered from all blocks.
// data: pairs {name, values}
// path: output path
template <class M>
void DumpCsv(
    const std::vector<std::pair<std::string, std::vector<typename M::Scal>>>&
        data,
    std::string path, M& m, char delim = ',');

} // namespace dump
