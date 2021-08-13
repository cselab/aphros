// Created by Petr Karnakov on 09.08.2021
// Copyright 2021 ETH Zurich

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>

#include "distr/distrbasic.h"
#include "dump/dump.h"
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
  struct {
    std::vector<Vect> x;
    std::vector<Scal> owner_init; // block owning particles after initialization
    std::vector<Scal> owner_comm; // block owning particles after communication
    std::vector<std::pair<std::string, std::vector<Scal>>> csvdata;
  } * ctx(sem);
  auto& t = *ctx;
  if (sem()) {
    std::mt19937 gen(17 + m.GetId());
    std::uniform_real_distribution<double> uniform;
    const int npart = var.Int["npart"];
    // Seed particles
    for (int i = 0; i < npart; ++i) {
      Vect x;
      for (auto d : m.dirs) {
        x[d] = uniform(gen);
      }
      t.x.push_back(x);
    }
    t.owner_init.resize(t.x.size(), m.GetId());
  }
  if (sem()) {
    typename M::CommPartRequest req;
    req.x = &t.x;
    req.attr_scal = {&t.owner_init};
    req.attr_vect = {};
    m.CommPart(req);
  }
  if (sem()) {
    t.owner_comm.resize(t.x.size(), m.GetId());
    // Add coordinates
    for (auto d : m.dirs) {
      const std::string name(1, M::Dir(d).letter());
      t.csvdata.emplace_back(name, std::vector<Scal>());
      for (auto x : t.x) {
        t.csvdata.back().second.push_back(x[d]);
      }
    }
    // Add other fields
    t.csvdata.emplace_back("owner_init", t.owner_init);
    t.csvdata.emplace_back("owner_comm", t.owner_comm);
  }
  if (sem.Nested()) {
    DumpCsv(t.csvdata, "particles.csv", m);
  }
  if (sem()) {
    DumpCsv(t.csvdata, util::Format("particles_{}.csv", m.GetId()));
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
  parser.AddVariable<int>("--npart", 32).Help("Number of particles per block");
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
  conf += "\nset int npart " + args.Int.GetStr("npart");
  conf += "\n" + args.String["extra"];

  return RunMpiBasicString<M>(mpi, Run, conf);
}
