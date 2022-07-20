// Created by Petr Karnakov on 26.01.2021
// Copyright 2021 ETH Zurich

#include <string>
#include <vector>

#include "dump/xmf.h"
#include "geom/vect.h"
#include "parse/argparse.h"
#include "util/format.h"

template <size_t dim_>
struct Mesh {
  static constexpr size_t dim = dim_;
  using Scal = double;
  using Vect = generic::Vect<Scal, dim>;
  using MIdx = generic::MIdx<dim>;
};

template <int dim_>
int Run(const Vars& args) {
  using M = Mesh<dim_>;
  constexpr size_t dim = M::dim;
  using MIdx = typename M::MIdx;
  using Vect = typename M::Vect;

  auto args_vect = [&](std::string key) {
    auto v = args.Vect[key];
    fassert(v.size() >= dim, ", incorrect size of vector '" + key + "'");
    return Vect(v.data());
  };
  auto args_midx = [&](std::string key) {
    auto v = args.Vect[key];
    fassert(v.size() >= dim, ", incorrect size of vector '" + key + "'");
    return MIdx(Vect(v.data()) + Vect(0.5));
  };

  using Xmf = dump::Xmf<Vect>;
  typename Xmf::Meta meta;
  meta.name = args.String["name"];
  meta.binpath = args.String["binpath"];
  meta.type = dump::StringToType(args.String["type"]);
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

  Xmf::WriteXmf(args.String["output"], meta, false);
  return 0;
}

int main(int argc, const char** argv) {
  ArgumentParser parser("Writes metadata");
  parser.AddVariable<std::string>("--name", "u").Help("Field name");
  parser.AddVariable<std::string>("--binpath", "u.raw")
      .Help("Path to binary datafile");
  parser.AddVariable<std::string>("--type", "Double")
      .Help("Number type")
      .Options({"UShort", "Float", "Double"});
  parser.AddVariable<int>("--dim", 3)
      .Help("Space dimensionality")
      .Options({1, 2, 3, 4});
  parser.AddVariable<std::vector<double>>("--dimensions", {8, 16, 24, 8})
      .Help("Full size, number of cells");
  parser.AddVariable<std::vector<double>>("--start", {0, 0, 0, 0})
      .Help("Hyperslab start");
  parser.AddVariable<std::vector<double>>("--stride", {1, 1, 1, 1})
      .Help("Hyperslab stride");
  parser.AddVariable<std::vector<double>>("--count", {0, 0, 0, 0})
      .Help("Hyperslab count. If zero, computed from dimensions");
  parser.AddVariable<int>("--seek", 0).Help("Number of bytes to skip");
  parser.AddVariable<std::vector<double>>("--origin", {0, 0, 0, 0})
      .Help("Position of the grid node with zero index");
  parser.AddVariable<std::vector<double>>("--spacing", {1, 1, 1, 1})
      .Help("Spacing between grid nodes");

  parser.AddVariable<std::string>("output", "u.xmf")
      .Help("Path to output metadata file");
  auto args = parser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  switch (args.Int["dim"]) {
    case 1:
      return Run<1>(args);
    case 2:
      return Run<2>(args);
    case 3:
      return Run<3>(args);
    case 4:
      return Run<4>(args);
    default:
      return 1;
  }
}
