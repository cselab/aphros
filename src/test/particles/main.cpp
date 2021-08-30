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

void Run(M& m, Vars& var) {
  auto sem = m.GetSem("Run");
  struct {
    std::vector<Vect> x;
    std::vector<Scal> x_first; // first component of position
    std::vector<Vect> x_copy;
    std::vector<bool> is_inner; // true if particle is owned by current block
    std::vector<std::string> msg;
    std::ofstream fout;
  } * ctx(sem);
  auto& t = *ctx;

  auto stat_to_string = [&]() -> std::string {
    std::string res;
    Vect x_mean(0);
    Vect x_copy_mean(0);
    Scal x_first_mean = 0;
    size_t npart_inner = 0;
    for (size_t i = 0; i < t.x.size(); ++i) {
      if (t.is_inner[i]) {
        x_mean += t.x[i];
        x_copy_mean += t.x_copy[i];
        x_first_mean += t.x_first[i];
        ++npart_inner;
      }
    }
    x_mean /= npart_inner;
    x_copy_mean /= npart_inner;
    x_first_mean /= npart_inner;
    res += util::Format(
        "block {:}, particles {:}, inner particles {:}, x_mean {:}, "
        "x_copy_mean {:}, x_first_mean {:}\n",
        m.GetId(), t.x.size(), npart_inner, x_mean, x_copy_mean, x_first_mean);
    return res;
  };

  if (sem()) {
    std::mt19937 gen(17 + m.GetId());
    std::uniform_real_distribution<double> uniform;
    const int npart = var.Int["npart"];
    // Seed particles uniformly throughout the domain,
    // even outside the current block.
    for (int i = 0; i < npart; ++i) {
      Vect x;
      for (auto d : m.dirs) {
        x[d] = uniform(gen);
      }
      t.x.push_back(x);
    }
    t.x_first.resize(t.x.size());
    for (size_t i = 0; i < t.x.size(); ++i) {
      t.x_first[i] = t.x[i][0];
    }
    t.x_copy = t.x;
    t.is_inner.resize(t.x.size(), true);

    if (m.IsRoot()) {
      t.fout.open("stat.log");
    }
  }
  if (sem()) {
    t.msg = {stat_to_string()};
    m.Reduce(&t.msg, Reduction::concat);

    typename M::CommPartRequest req;
    req.x = &t.x;
    req.is_inner = &t.is_inner;
    req.attr_scal = {&t.x_first};
    req.attr_vect = {&t.x_copy};
    m.CommPart(req);
  }
  if (sem()) {
    if (m.IsRoot()) {
      std::sort(t.msg.begin(), t.msg.end());
      t.fout << "Before communication:\n";
      for (auto& s : t.msg) {
        t.fout << s;
      }
    }
  }
  if (sem()) {
    t.msg = {stat_to_string()};
    m.Reduce(&t.msg, Reduction::concat);
  }
  if (sem()) {
    if (m.IsRoot()) {
      std::sort(t.msg.begin(), t.msg.end());
      t.fout << "After communication:\n";
      for (auto& s : t.msg) {
        t.fout << s;
      }
    }
  }
}

int main(int argc, const char** argv) {
  MpiWrapper mpi(&argc, &argv);

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
  parser.AddVariable<int>("--npart", 128).Help("Number of particles per block");
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
