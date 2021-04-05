// Created by Petr Karnakov on 25.09.2020
// Copyright 2020 ETH Zurich

#include <cassert>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include <distr/distrbasic.h>
#include <dump/hdf.h>
#include <func/init.h>
#include <parse/argparse.h>
#include <parse/vars.h>
#include <util/distr.h>
#include <util/filesystem.h>
#include <util/format.h>
#include <util/vof.h>

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

void Init(
    GRange<size_t> layers, Multi<FieldCell<Scal>>& fcu,
    std::string prefix, const Vars& var, M& m) {
  auto sem = m.GetSem();
  struct {
    Vars var;
  } * ctx(sem);
  auto& t = *ctx;
  if (sem()) {
    fcu.Reinit(layers, m, 0);
  }
  for (auto l : layers) {
    const auto sl = std::to_string(l);
    if (sem()) {
      t.var.String.Set("init_vf", "list");
      t.var.String.Set("list_path", var.String["init_" + prefix + sl]);
      t.var.Int.Set("dim", m.GetGlobalSize()[2] == 1 ? 2 : 3);
      t.var.Int.Set("list_ls", 1);
    }
    if (sem.Nested()) {
      InitVf(fcu[l], t.var, m, false);
    }
  }
}

void Dump(
    GRange<size_t> layers, const Multi<FieldCell<Scal>>& fcu,
    std::string prefix, M& m) {
  auto sem = m.GetSem();
  for (auto l : layers) {
    const auto sl = std::to_string(l);
    const auto name = prefix + sl;
    const auto path = name + ".h5";
    if (sem.Nested()) {
      Hdf<M>::Write(fcu[l], path, m);
    }
    if (sem() && m.IsRoot()) {
      Hdf<M>::WriteXmf(util::SplitExt(path)[0] + ".xmf", prefix, path, m);
    }
  }
}

void Run(M& m, Vars& var) {
  auto sem = m.GetSem();
  struct {
    GRange<size_t> layers;
    Multi<FieldCell<Scal>> fcvf; // volume fraction
    Multi<FieldCell<Scal>> fccl; // color
    MapEmbed<BCond<Scal>> mebc; // boundary conditions for volume fraction
  } * ctx(sem);
  constexpr Scal kClNone = -1;
  auto& t = *ctx;
  if (sem()) {
    t.layers = GRange<size_t>(var.Int["layers"]);
    t.fccl.Reinit(t.layers, m, 0);
  }
  if (sem.Nested()) {
    Init(t.layers, t.fcvf, "vf", var, m);
  }
  if (sem()) {
    // Clear color in cells with zero volume fraction
    for (auto l : t.layers) {
      for (auto c : m.AllCells()) {
        t.fccl[l][c] = (t.fcvf[l][c] > 0 ? 1 : kClNone);
      }
    }
  }
  if (sem.Nested()) {
    const bool verbose = false;
    const bool reduce = false;
    const bool unionfind = false;
    const bool grid = false;
    UVof<M>().Recolor(
        t.layers, t.fcvf, t.fccl, t.fccl, 0, Vect(0), 1e10, t.mebc, verbose,
        unionfind, reduce, grid, m);
  }
  if (sem.Nested()) {
    Dump(t.layers, t.fcvf, "vf", m);
  }
  if (sem.Nested()) {
    Dump(t.layers, t.fccl, "cl", m);
  }
}

int main(int argc, const char** argv) {
  MpiWrapper mpi(&argc, &argv);

  ArgumentParser parser(
      "Example for connected component labeling", mpi.IsRoot());
  parser.AddVariable<std::string>("--extra", "")
      .Help("Extra configuration (commands 'set ... ')");
  parser.AddVariable<std::string>("config", "a.conf")
      .Help("Path to configuration file");
  parser.AddVariable<int>("--nx", 128).Help("Mesh size in the x-direction");
  parser.AddVariable<int>("--ny", 128).Help("Mesh size in the y-direction");
  parser.AddVariable<int>("--nz", 1).Help("Mesh size in the z-direction");
  parser.AddVariable<int>("--bs", 32).Help("Block size");
  parser.AddVariable<int>("--layers", 4).Help("Number of layers");
  auto args = parser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  std::stringstream conf;

  MIdx meshsize(args.Int["nx"], args.Int["ny"], args.Int["nz"]);
  const int bs = args.Int["bs"];
  MIdx blocksize(bs, bs, args.Int["nz"] == 1 ? 1 : bs);
  Subdomains<MIdx> sub(meshsize, blocksize, mpi.GetCommSize());
  conf << sub.GetConfig() << '\n';

  const auto configpath = args.String["config"];
  if (!configpath.empty()) {
    conf << "include " + configpath + '\n';
  }
  conf << "set string backend native\n";
  conf << "set double extent 1\n";
  conf << util::Format("set int layers {}\n", args.Int["layers"]);
  conf << args.String["extra"] << '\n';

  return RunMpiBasicString<M>(mpi, Run, conf.str());
}
