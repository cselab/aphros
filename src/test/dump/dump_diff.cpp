// Created by Petr Karnakov on 26.01.2021
// Copyright 2021 ETH Zurich

#include <omp.h>
#include <iostream>

#include "debug/linear.h"
#include "distr/distrbasic.h"
#include "dump/hdf.h"
#include "dump/raw.h"
#include "parse/argparse.h"
#include "util/distr.h"
#include "util/filesystem.h"
#include "util/format.h"

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

std::string GetFormat(std::string path, std::string format) {
  const auto ext = util::SplitExt(path)[1];
  if (format == "auto") {
    if (ext == ".h5") {
      format = "h5";
    } else if (ext == ".raw") {
      format = "raw";
    } else {
      fassert(false, util::Format("Unknown extension '{}' of '{}'", ext, path));
    }
  }
  return format;
}

void Run(M& m, Vars& var) {
  auto sem = m.GetSem(__func__);
  using Raw = dump::Raw<M>;
  struct {
    FieldCell<Scal> fc_read1;
    FieldCell<Scal> fc_read2;
    FieldCell<Scal> fc_diff;
    std::vector<generic::Vect<Scal, 3>> norms;
    Raw::Meta meta;
  } * ctx(sem);
  auto& t = *ctx;

  auto read = [&](FieldCell<Scal>& fc_buf, std::string path) {
    if (sem("readxmf")) {
      auto format = GetFormat(path, var.String["format"]);
      if (format == "raw") {
        const auto xmfpath = util::SplitExt(path)[0] + ".xmf";
        t.meta = dump::Xmf<Vect>::ReadXmf(xmfpath);
      }
    }
    if (sem.Nested("read")) {
      auto format = GetFormat(path, var.String["format"]);
      if (format == "h5") {
        dump::Hdf<M>::Read(fc_buf, path, m);
      } else if (format == "raw") {
        Raw::Read(fc_buf, t.meta, path, m);
      } else {
        fassert(false, "Unkown format=" + format);
      }
    }
  };

  if (sem()) {
    t.fc_read1.Reinit(m, 0);
    t.fc_read2.Reinit(m, 0);
  }
  read(t.fc_read1, var.String["input1"]);
  read(t.fc_read2, var.String["input2"]);
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

  ArgumentParser parser(
      "Compares two scalar arrays in files"
      ". Prints the L1,L2, and Linf norms of the difference"
      ". Assumes both arrays have the same shape",
      mpi.IsRoot());
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

  auto get_shape = [&](std::string path) {
    MIdx res;
    const auto format = GetFormat(path, args.String.GetStr("format"));
    if (format == "h5") {
      const auto shape = dump::Hdf<M>::GetShape(path);
      res[0] = shape[2];
      res[1] = shape[1];
      res[2] = shape[0];
    } else if (format == "raw") {
      const auto xmfpath = util::SplitExt(path)[0] + ".xmf";
      const auto meta = dump::Xmf<Vect>::ReadXmf(xmfpath);
      res = meta.count;
    } else {
      fassert(
          false, "Can't determine the field dimensions of format '" + format +
                     "' of '" + path + "'");
    }
    return res;
  };

  const MIdx mesh_size = get_shape(args.String.GetStr("input1"));
  const MIdx mesh_size2 = get_shape(args.String.GetStr("input2"));
  fassert_equal(mesh_size, mesh_size2, ". Files have different shapes");

  MIdx block_size(args.Int["bs"]);
  if (mesh_size[2] == 1) {
    block_size[2] = 1;
  }

  fassert(
      mesh_size % block_size == MIdx(0),
      util::Format(
          "File shape {} not divisible by block size {}", mesh_size,
          block_size));

  Subdomains<MIdx> sub(mesh_size, block_size, mpi.GetCommSize());
  conf += sub.GetConfig();
  conf += "\nset string input1 " + args.String.GetStr("input1");
  conf += "\nset string input2 " + args.String.GetStr("input2");
  conf += "\nset string format " + args.String.GetStr("format");

  return RunMpiBasicString<M>(mpi, Run, conf);
}
