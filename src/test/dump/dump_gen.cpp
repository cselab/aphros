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

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

void Run(M& m, Vars& var) {
  auto sem = m.GetSem(__func__);
  struct {
    FieldCell<Scal> fc_write;
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
    const auto content = var.String["content"];
    t.fc_write.Reinit(m, 0);
    if (content == "norm") {
      for (auto c : m.CellsM()) {
        t.fc_write[c] = c.center().norm();
      }
    } else if (content == "zero") {
      // nop
    } else if (content == "one") {
      t.fc_write.Reinit(m, 1);
    } else {
      fassert(false, "Unknown content=" + content);
    }
  }
  const auto output = var.String["output_path"];
  auto format = get_format(output);
  if (format == "h5") {
    if (sem.Nested("write")) {
      Hdf<M>::Write(t.fc_write, output, m);
    }
    if (sem("writexmf")) {
      Hdf<M>::WriteXmf(util::SplitExt(output)[0] + ".xmf", "u", output, m);
    }
  } else if (format == "raw") {
    dump::Raw<M>::Meta meta;
    meta.size = m.GetGlobalSize();
    meta.count = m.GetGlobalSize();
    using Raw = dump::Raw<M>;
    if (sem.Nested("write")) {
      Raw::Write(t.fc_write, meta, output, m);
    }
    if (sem("writexmf")) {
      Raw::WriteXmf(
          util::SplitExt(output)[0] + ".xmf", "u", Raw::Type::Float64, output,
          m);
    }
  } else {
    fassert(false, "Unkown format=" + format);
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
  parser.AddVariable<std::string>({"--content", "-c"}, "norm")
      .Help("Field to write")
      .Options({"norm", "zero", "one"});
  parser.AddVariable<std::string>("output", "o.h5").Help("Path to output file");
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
  // XXX Using "output_path" since "output" is aleady used in default config
  conf += "\nset string output_path " + args.String.GetStr("output");
  conf += "\nset string format " + args.String.GetStr("format");
  conf += "\nset string content " + args.String.GetStr("content");

  return RunMpiBasicString<M>(mpi, Run, conf);
}
