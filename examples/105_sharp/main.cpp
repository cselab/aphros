// Created by Petr Karnakov on 25.09.2020
// Copyright 2020 ETH Zurich

#include <cassert>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include <distr/distrbasic.h>
#include "dump/hdf.h"
#include "dump/raw.h"
#include "parse/argparse.h"
#include "parse/vars.h"
#include "parse/vof.h"
#include "solver/vof.h"
#include "util/distr.h"
#include "util/filesystem.h"
#include "util/hydro.h"
#include "util/vof.h"

using M = MeshStructured<double, 3>;
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
  const bool verbose = var.Int["VERBOSE"];
  using Scal = typename M::Scal;
  auto sem = m.GetSem();
  struct {
    std::unique_ptr<Vof<M>> solver;
    Vof<M>::Par par;

    FieldCell<Scal> fcu;
    FieldEmbed<Scal> fe_flux;
    FieldCell<Scal> fc_src;
    FieldCell<Scal> fc_zero;
    FieldCell<Vect> fc_zerov;

    // boundary conditions for advection (empty)
    MapEmbed<BCondAdvection<Scal>> bc;

    // for csv dump
    UVof<M> uvof;
    FieldCell<Scal> fccl;
    std::vector<std::string> column_names;
    std::vector<Scal> row_colors;
    std::vector<std::vector<Scal>> table;
  } * ctx(sem);
  auto& t = *ctx;

  auto read = [&](FieldCell<Scal>& fc_buf, std::string path) {
    if (sem.Nested("read")) {
      auto format = GetFormat(path, var.String["format"]);
      if (format == "h5") {
        Hdf<M>::Read(fc_buf, path, m);
      } else if (format == "raw") {
        dump::Raw<M>::Meta meta;
        meta.size = m.GetGlobalSize();
        meta.count = m.GetGlobalSize();
        using Raw = dump::Raw<M>;
        Raw::Read(fc_buf, meta, path, m);
      } else {
        fassert(false, "Unkown format=" + format);
      }
    }
  };
  auto write = [&](auto fc_buf, std::string path) {
    auto format = GetFormat(path, var.String["format"]);
    if (format == "h5") {
      if (sem.Nested("write")) {
        Hdf<M>::Write(fc_buf(), path, m);
      }
      if (sem("writexmf")) {
        Hdf<M>::WriteXmf(util::SplitExt(path)[0] + ".xmf", "u", path, m);
      }
    } else if (format == "raw") {
      dump::Raw<M>::Meta meta;
      meta.size = m.GetGlobalSize();
      meta.count = m.GetGlobalSize();
      using Raw = dump::Raw<M>;
      if (sem.Nested("write")) {
        Raw::Write(fc_buf(), meta, path, m);
      }
      if (sem("writexmf")) {
        Raw::WriteXmf(
            util::SplitExt(path)[0] + ".xmf", "u", Raw::Type::Float64, path, m);
      }
    } else {
      fassert(false, "Unkown format=" + format);
    }
  };
  read(t.fcu, var.String["inputpath"]);
  if (sem("ctor")) {
    t.fe_flux.Reinit(m, 0);
    t.fc_src.Reinit(m, 0);
    t.fc_zero.Reinit(m, 0);
    t.fc_zerov.Reinit(m, Vect(0));

    t.par.clipth = 1e-10;
    t.par.sharpen = true;
    t.par.sharpen_cfl = var.Double["cfl"];
    t.par.dim = var.Int["dim"];
    const FieldCell<Scal> fccl(m, 0);
    t.solver.reset(new Vof<M>(
        m, m, t.fcu, fccl, t.bc, &t.fe_flux, &t.fc_src, 0., 1., t.par));
  }
  if (sem.Nested("start")) {
    t.solver->StartStep();
  }
  for (int i = 0; i < var.Int["steps"]; ++i) {
    if (verbose && sem()) {
      if (m.IsRoot()) {
        std::cout << "step " << i << std::endl;
      }
    }
    if (i == 0) {
      if (sem.Nested("iter")) {
        t.solver->MakeIteration();
      }
    } else {
      if (sem.Nested("sharp")) {
        t.solver->Sharpen();
      }
    }
  }
  if (sem.Nested("finish")) {
    t.solver->FinishStep();
  }
  write(
      [&]() { return t.solver->GetField(Step::iter_curr); },
      var.String["outputpath"]);
  const std::string* vtk_out = var.String.Find("vtk_out");
  if (vtk_out && sem.Nested()) {
    t.solver->DumpInterface(*vtk_out);
  }
  const std::string* vtk_out_march = var.String.Find("vtk_out_march");
  if (vtk_out_march && sem.Nested()) {
    t.solver->DumpInterfaceMarch(*vtk_out_march);
  }
  const std::string* csv_out = var.String.Find("csv_out");
  if (csv_out && sem()) {
    t.fccl.Reinit(m, 0);
    for (auto c : m.SuCells()) {
      t.fccl[c] = t.solver->GetField()[c] > 0 ? 0 : Vof<M>::kClNone;
    }
  }
  if (csv_out && sem.Nested()) {
    auto plic = t.solver->GetPlic();
    t.uvof.Recolor(
        plic.layers, plic.vfcu, &t.fccl, &t.fccl, t.par.clfixed,
        t.par.clfixed_x, t.par.coalth, {}, t.par.verb, t.par.recolor_unionfind,
        t.par.recolor_reduce, t.par.recolor_grid, m);
  }
  if (csv_out && sem.Nested()) {
    auto plic = t.solver->GetPlic();
    CalcTraj<M>(
        m, plic.layers, plic.vfcu, &t.fccl, plic.vfcim, t.fc_zero, t.fc_zerov,
        t.column_names, t.row_colors, t.table);
  }
  if (csv_out && sem()) {
    if (m.IsRoot()) {
      std::ofstream o;
      o.open(*csv_out);
      o.precision(16);
      // header
      {
        o << "cl";
        for (size_t i = 0; i < t.column_names.size(); ++i) {
          o << "," << t.column_names[i];
        }
        o << std::endl;
      }
      // content
      for (size_t i = 0; i < t.row_colors.size(); ++i) {
        o << t.row_colors[i];
        for (auto v : t.table[i]) {
          o << "," << v;
        }
        o << "\n";
      }
    }
  }
  if (sem()) { // XXX empty stage
  }
}

int main(int argc, const char** argv) {
  MpiWrapper mpi(&argc, &argv);

  ArgumentParser parser(
      "Sharpens the image using PLIC advection", mpi.IsRoot());
  parser.AddSwitch({"--verbose", "-v"}).Help("Report steps");
  parser.AddVariable<double>("--cfl", 0.1)
      .Help("CFL number for advection, valid values between 0 and 1");
  parser.AddVariable<std::string>("--vtk_out_march")
      .Help("Path to output VTK with surface from marching cubes");
  parser.AddVariable<std::string>("--vtk_out")
      .Help("Path to output VTK with piecewise linear surface");
  parser.AddVariable<std::string>("--csv_out")
      .Help("Path to output CSV with centroids of connected components");
  parser.AddVariable<int>("--steps", 5).Help("Number of sharpening steps");
  parser.AddVariable<std::string>("--extra", "")
      .Help("Extra configuration (commands 'set ... ')");
  parser.AddVariable<std::string>({"--format", "-f"}, "auto")
      .Help("File format")
      .Options({"auto", "h5", "raw", "dat"});

  parser.AddVariable<std::string>("input").Help(
      "Path to input image as HDF5 array of floats between 0 and 1");
  parser.AddVariable<std::string>("output").Help(
      "Path to output image, can be the same as input");

  auto args = parser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  int nx, ny, nz;
  const auto inputformat =
      GetFormat(args.String.GetStr("input"), args.String.GetStr("format"));
  if (inputformat == "h5") {
    const auto shape = Hdf<M>::GetShape(args.String["input"]);
    nx = shape[2];
    ny = shape[1];
    nz = shape[0];
  } else {
    nx = 0;
    ny = 0;
    nz = 0;
  }

  auto checkdiv = [](int n, int b, std::string name) {
    if (n % b) {
      std::stringstream s;
      s << name << "=" << n << " not divisible by " << b;
      throw std::runtime_error(s.str());
    }
  };

  int bs = 8;
  int bsz = 8;

  if (nz == 1) { // 2D
    checkdiv(nx, bs, "NX");
    checkdiv(ny, bs, "NY");
    bsz = 1;
    while (nx % bs == 0 && ny % bs == 0 && bs <= 32) {
      bs *= 2;
    }
    bs /= 2;
  } else { // 3D
    checkdiv(nx, bs, "NX");
    checkdiv(ny, bs, "NY");
    checkdiv(nz, bsz, "NZ");
    while (nx % bs == 0 && ny % bs == 0 && nz % bs == 0 && bs <= 32) {
      bs *= 2;
    }
    bs /= 2;
    bsz = bs;
  }
  int bsy = bs;

  std::stringstream conf;
  Subdomains<MIdx> sub(MIdx(nx, ny, nz), MIdx(bs, bsy, bsz), mpi.GetCommSize());
  conf << sub.GetConfig() << '\n';

  conf << "set int dim " << (nz == 1 ? 2 : 3) << '\n';

  conf << "set int steps " << args.Int["steps"] << '\n';

  conf << "set string inputpath " << args.String["input"] << '\n';
  conf << "set string outputpath " << args.String["output"] << '\n';
  if (auto* p = args.String.Find("vtk_out")) {
    conf << "set string vtk_out " << *p << '\n';
  }
  if (auto* p = args.String.Find("vtk_out_march")) {
    conf << "set string vtk_out_march " << *p << '\n';
  }
  if (auto* p = args.String.Find("csv_out")) {
    conf << "set string csv_out " << *p << '\n';
  }
  conf << "set double cfl " << args.Double["cfl"] << '\n';
  conf << "set int VERBOSE " << args.Int["verbose"] << '\n';
  conf << "set string format " << args.String["format"] << '\n';

  conf << args.String["extra"] << '\n';

  return RunMpiBasicString<M>(mpi, Run, conf.str());
}
