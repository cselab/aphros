// Created by Petr Karnakov on 26.01.2021
// Copyright 2021 ETH Zurich

#include <omp.h>
#include <iostream>

#include "debug/linear.h"
#include "distr/distrbasic.h"
#include "dump/hdf.h"
#include "parse/argparse.h"
#include "util/distr.h"
#include "util/filesystem.h"
#include "util/format.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

void Run(M& m, Vars& var) {
  auto sem = m.GetSem(__func__);
  struct {
    FieldCell<Scal> fc_read1;
    FieldCell<Scal> fc_read2;
    FieldCell<Scal> fc_diff;
    std::vector<generic::Vect<Scal, 3>> norms;
  } * ctx(sem);
  auto& t = *ctx;
  auto get_format = [&var](std::string path) {
    const auto ext = util::SplitExt(path)[1];
    auto format = var.String["format"];
    if (format == "auto") {
      if (ext == ".h5") {
        format = "h5";
      } else if (ext == ".raw") {
        format = "raw";
      } else if (ext == ".dat") {
        format = "dat";
      } else {
        fassert(
            false, util::Format("Unknown extension '{}' of '{}'", ext, path));
      }
    }
    return format;
  };
  if (sem()) {
    t.fc_read1.Reinit(m, 0);
    t.fc_read2.Reinit(m, 0);
  }
  auto read = [&](FieldCell<Scal>& fc_buf, std::string suff) {
    if (sem.Nested("read" + suff)) {
      const auto input = var.String["input" + suff];
      auto format = get_format(input);
      if (format == "h5") {
        Hdf<M>::Read(fc_buf, input, m);
      } else {
        fassert(false, "Unkown format=" + format);
      }
    }
  };
  read(t.fc_read1, "1");
  read(t.fc_read2, "2");
  if (sem("norms")) {
    t.fc_diff.Reinit(m, 0);
    for (auto c : m.Cells()) {
      t.fc_diff[c] = t.fc_read2[c] - t.fc_read1[c];
    }
  }
  if (sem.Nested("norms")) {
    t.norms = UDebug<M>::GetNorms({&t.fc_diff}, m);
  }
  if (sem("print")) {
    if (m.IsRoot()) {
      const auto& n = t.norms[0];
      std::cout << n[0] << ' ' << n[1] << ' ' << n[2] << '\n';
    }
  }
  if (sem()) {
  }
}

int main(int argc, const char** argv) {
  MpiWrapper mpi(&argc, &argv);

  ArgumentParser parser("Test for writers and readers.", mpi.IsRoot());
  parser.AddVariable<int>("--nx", 16).Help("Mesh size in x-direction");
  parser.AddVariable<int>("--ny", 16).Help("Mesh size in y-direction");
  parser.AddVariable<int>("--nz", 16).Help("Mesh size in z-direction");
  parser.AddVariable<int>("--bs", 8)
      .Help("Block size in all directions")
      .Options({8, 16, 32});
  parser.AddVariable<std::string>({"--format", "-f"}, "auto")
      .Help("File format")
      .Options({"auto", "h5", "raw", "dat"});
  parser.AddVariable<std::string>("input1").Help("Path to input file 1");
  parser.AddVariable<std::string>("input2").Help("Path to input file 2");
  auto args = parser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  std::string conf;

  MIdx mesh_size(args.Int["nx"], args.Int["ny"], args.Int["nz"]);
  MIdx block_size(args.Int["bs"]);
  if (mesh_size[2] == 1) {
    block_size[2] = 1;
  }

  Subdomains<MIdx> sub(mesh_size, block_size, mpi.GetCommSize());
  conf += sub.GetConfig();
  conf += "\nset string input1 " + args.String.GetStr("input1");
  conf += "\nset string input2 " + args.String.GetStr("input2");
  conf += "\nset string format " + args.String.GetStr("format");

  return RunMpiBasicString<M>(mpi, Run, conf);
}
