// Created by Petr Karnakov on 26.01.2021
// Copyright 2021 ETH Zurich

#include <string>
#include <vector>

#include "dump/raw.h"
#include "dump/raw.ipp"
#include "geom/vect.h"
#include "parse/argparse.h"
#include "util/format.h"

struct M {
  static constexpr size_t dim = 3;
  using Scal = double;
  using Vect = generic::Vect<Scal, dim>;
  using MIdx = generic::Vect<int, dim>;
};

using MIdx = M::MIdx;
using Vect = M::Vect;
using Scal = M::Scal;

int main(int argc, const char** argv) {
  MpiWrapper mpi(&argc, &argv);

  ArgumentParser parser("Writes metadata.", mpi.IsRoot());
  parser.AddVariable<std::string>("--name", "u").Help("Field name");
  parser.AddVariable<std::string>("--binpath", "u.raw")
      .Help("Path to binary datafile");
  parser.AddVariable<std::string>("--type", "Double")
      .Help("Number type")
      .Options({"UShort", "Float", "Double"});
  parser.AddVariable<std::vector<double>>("--dimensions", {8, 16, 24})
      .Help("Full size, number of cells");
  parser.AddVariable<std::vector<double>>("--start", {0, 0, 0})
      .Help("Hyperslab start");
  parser.AddVariable<std::vector<double>>("--stride", {1, 1, 1})
      .Help("Hyperslab stride");
  parser.AddVariable<std::vector<double>>("--count", {0, 0, 0})
      .Help("Hyperslab count. If zero, computed from dimensions");
  parser.AddVariable<int>("--seek", 0).Help("Number of bytes to skip");
  parser.AddVariable<std::vector<double>>("--origin", {0, 0, 0})
      .Help("Position of the grid node with zero index");
  parser.AddVariable<std::vector<double>>("--spacing", {1, 1, 1})
      .Help("Spacing between grid nodes");

  parser.AddVariable<std::string>("output", "u.xmf")
      .Help("Path to output metadata file");
  auto args = parser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  auto args_vect = [&](std::string key) {
    auto v = args.Vect[key];
    fassert_equal(v.size(), 3, ", incorrect size of vector '" + key + "'");
    return Vect(v.data());
  };
  auto args_midx = [&](std::string key) {
    auto v = args.Vect[key];
    fassert_equal(v.size(), 3, ", incorrect size of vector '" + key + "'");
    return MIdx(Vect(v.data()) + Vect(0.5));
  };

  using Raw = dump::Raw<M>;
  Raw::Meta meta;
  meta.name = args.String["name"];
  meta.binpath = args.String["binpath"];
  meta.type = Raw::StringToType(args.String["type"]);
  meta.dimensions = args_midx("dimensions");
  meta.start = args_midx("start");
  meta.stride = args_midx("stride");
  meta.count = args_midx("count");
  if (meta.count == MIdx(0)) {
    meta.count = (meta.dimensions - meta.start) / meta.stride;
  }
  meta.seek = args.Int["seek"];
  meta.origin = args_vect("origin");
  meta.spacing = args_vect("spacing");

  Raw::WriteXmf("o.xmf", meta);
}
