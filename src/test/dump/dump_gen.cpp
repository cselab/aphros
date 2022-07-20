// Created by Petr Karnakov on 26.01.2021
// Copyright 2021 ETH Zurich

#include <omp.h>
#include <iostream>

#include "debug/linear.h"
#include "distr/distrbasic.h"
#include "dump/hdf.h"
#include "dump/raw.h"
#include "dump/xmf.h"
#include "parse/argparse.h"
#include "util/distr.h"
#include "util/filesystem.h"
#include "util/format.h"

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

void Run(M& m, Vars& var) {
  auto sem = m.GetSem(__func__);
  using Raw = dump::Raw<M>;
  using Xmf = dump::Xmf<Vect>;
  struct {
    FieldCell<Scal> fc_write;
    Xmf::Meta meta;
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
      dump::Hdf<M>::Write(t.fc_write, output, m);
    }
    if (sem("writexmf") && m.IsRoot()) {
      dump::Hdf<M>::WriteXmf(
          util::SplitExt(output)[0] + ".xmf", "u", output, m);
    }
  } else if (format == "raw") {
    if (sem("writexmf")) {
      t.meta = Xmf::GetMeta(m);
      t.meta.name = "u";
      t.meta.binpath = output;
      const auto type = var.String["type"];
      if (type.length()) {
        t.meta.type = dump::StringToType(type);
      } else {
        t.meta.type = dump::Type::Float64;
      }
      if (m.IsRoot()) {
        Xmf::WriteXmf(util::SplitExt(output)[0] + ".xmf", t.meta, false);
      }
      if (t.meta.type == dump::Type::UInt16) {
        for (auto c : m.Cells()) {
          auto& u = t.fc_write[c];
          u = std::min(1., std::max(0., u)) *
              std::numeric_limits<std::uint16_t>::max();
        }
      }
    }
    if (sem.Nested("write")) {
      Raw::WriteMeshBlocks(t.fc_write, t.meta, output, m);
    }
  } else {
    fassert(false, "Unkown format=" + format);
  }
  if (sem()) {
  }
}

int main(int argc, const char** argv) {
  MpiWrapper mpi(&argc, &argv);

  ArgumentParser parser("Writes a file with scalar field.", mpi.IsRoot());
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
  parser.AddVariable<std::string>("--type", "")
      .Help("Number type of the output file. UShort is rescaled to [0, 1].")
      .Options({"", "UShort", "Float", "Double"});
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
  conf += "\nset string type " + args.String.GetStr("type");

  return RunMpiBasicString<M>(mpi, Run, conf);
}
