// Created by Petr Karnakov on 11.11.2020
// Copyright 2020 ETH Zurich

#ifdef _OPENMP
#include <omp.h>
#else
int omp_get_max_threads() { return 1; }
#endif

#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "distr/distrbasic.h"
#include "geom/mesh.h"
#include "kernel/kernelmeshpar.h"
#include "parse/argparse.h"
#include "util/distr.h"
#include "util/format.h"
#include "util/mpi.h"

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

static std::ofstream fout;

void Run(M& m, Vars& var) {
  auto sem = m.GetSem(__func__);
  if (sem()) {
    if (m.IsLead()) {
      fout.open(
          util::Format("out_{}", MpiWrapper(m.GetMpiComm()).GetCommRank()));
    }
  }
  if (sem()) {
    fout << util::Format(
        "GetMpiRankFromId() from mesh with id={:} rank={:}\n", m.GetId(),
        MpiWrapper(m.GetMpiComm()).GetCommRank());
    const MIdx procs(var.Int["px"], var.Int["py"], var.Int["pz"]);
    const MIdx blocks(var.Int["bx"], var.Int["by"], var.Int["bz"]);
    const int nblocks = (procs * blocks).prod();
    for (int id = 0; id < nblocks; ++id) {
      fout << util::Format("id={:} rank={:}\n", id, m.GetMpiRankFromId(id));
    }
  }
}

int main(int argc, const char** argv) {
  MpiWrapper mpi(&argc, &argv);

  ArgumentParser parser("Test for communication maps.", mpi.IsRoot());
  parser.AddVariable<int>("--nx", 16).Help("Mesh size in x-direction");
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
  const int block = args.Int["block"];
  Subdomains<MIdx> sub(
      MIdx(nx), MIdx(block), mpi.GetCommSize() / omp_get_max_threads());

  std::string conf = "";
  conf += sub.GetConfig();
  conf += "\n" + args.String["extra"];

  return RunMpiBasicString<M>(mpi, Run, conf);
}
