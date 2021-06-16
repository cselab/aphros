// Created by Petr Karnakov on 11.11.2020
// Copyright 2020 ETH Zurich

#include <omp.h>
#include <algorithm>
#include <iostream>
#include <map>

#include "distr/commmap.h"
#include "distr/distrbasic.h"
#include "parse/argparse.h"
#include "util/distr.h"
#include "util/timer.h"

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

void Run(M& m, Vars&) {
  auto sem = m.GetSem(__func__);
  using Expr = typename M::Expr;
  struct {
    CommMap<M> comm_map;
    FieldCell<Expr> fc_system;
  } * ctx(sem);

  auto& t = *ctx;

  if (sem()) {
    t.fc_system.Reinit(m);
    int i = 0;
    for (auto c : m.Cells()) {
      t.fc_system[c][0] = m.GetId() * m.GetInBlockCells().size() + (i++);
    }
  }
  if (sem.Nested()) {
    t.comm_map.Init(m);
  }
  if (sem.Nested()) {
    t.comm_map.ConvertSystem(t.fc_system, m);
  }
  if (sem.Nested()) {
    t.comm_map.PrintStat(m, std::cout);
  }
}

int main(int argc, const char** argv) {
  MpiWrapper mpi(&argc, &argv);

  ArgumentParser parser("Test for communication maps.", mpi.IsRoot());
  parser.AddVariable<int>("--nx", 16).Help("Mesh size in x-direction");
  parser.AddVariable<int>("--ny", 16).Help("Mesh size in y-direction");
  parser.AddVariable<int>("--block", 8)
      .Help("Block size in x- and y-directions")
      .Options({8, 16, 32});
  parser.AddVariable<std::string>("--extra", "")
      .Help("Extra configuration (commands 'set ... ')");
  auto args = parser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  const int nx = args.Int["nx"];
  const int ny = args.Int["ny"];
  const int block = args.Int["block"];
  Subdomains<MIdx> sub(
      MIdx(nx, ny, 1), MIdx(block, block, 1),
      mpi.GetCommSize() / omp_get_max_threads());

  std::string conf = "";
  conf += sub.GetConfig();
  conf += "\n" + args.String["extra"];

  return RunMpiBasicString<M>(mpi, Run, conf);
}
