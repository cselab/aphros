// Created by Petr Karnakov on 09.08.2021
// Copyright 2021 ETH Zurich

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#include "distr/distrbasic.h"
#include "geom/mesh.h"
#include "kernel/kernelmeshpar.h"
#include "parse/argparse.h"
#include "solver/pois.h"
#include "solver/solver.h"
#include "util/format.h"
#include "util/mpi.h"
#include "util/suspender.h"

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

template <class T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) {
  std::string p = "";
  for (auto a : v) {
    o << p << a;
    p = " ";
  }
  return o;
}

std::ofstream fout;

void Run(M& m, Vars& var) {
  auto sem = m.GetSem("Run");
  if (sem()) {
    if (m.IsLead()) {
    }
  }
  if (sem()) {
    std::vector<Vect> x;
    std::vector<Scal> sa;
    std::vector<Vect> va;
    typename M::CommPartRequest req;
    req.x = &x;
    req.attr_scal = {&sa};
    req.attr_vect = {&va};
    m.CommPart(req);
  }
}

int main(int argc, const char** argv) {
  MpiWrapper mpi(&argc, &argv);
  fout.open(util::Format("out_{}", mpi.GetCommRank()));

  ArgumentParser parser("Test for particle communication", mpi.IsRoot());
  auto instances = []() {
    auto map = ModuleDistr<M>::GetInstances();
    std::vector<std::string> res;
    for (auto p : map) {
      res.push_back(p.first);
    }
    return res;
  }();
  parser.AddVariable<std::string>("--backend", "native")
      .Help("Communication backend")
      .Options(instances);
  parser.AddVariable<int>("--block", 16).Help("Block size in all directions");
  parser.AddVariable<int>("--mesh", 32).Help("Mesh size in all directions");
  parser.AddVariable<std::string>("--extra", "")
      .Help("Extra configuration (commands 'set ... ')");
  auto args = parser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  std::string conf;

  MIdx mesh_size(args.Int["mesh"]);
  MIdx block_size(args.Int["block"]);

  Subdomains<MIdx> sub(mesh_size, block_size, mpi.GetCommSize());
  conf += sub.GetConfig();

  conf += "\nset string backend " + args.String["backend"];
  conf += "\n" + args.String["extra"];

  return RunMpiBasicString<M>(mpi, Run, conf);
}
