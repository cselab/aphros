// Created by Petr Karnakov on 25.09.2020
// Copyright 2020 ETH Zurich

#include <mpi.h>
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
#include "parse/argparse.h"
#include "parse/vars.h"
#include "parse/vof.h"
#include "solver/vof.h"
#include "util/hydro.h"
#include "util/vof.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

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
  if (sem.Nested()) {
    Hdf<M>::Read(t.fcu, var.String["hdf_in"], m);
  }
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
  if (sem.Nested()) {
    Hdf<M>::Write(
        t.solver->GetField(Step::iter_curr), var.String["hdf_out"], m);
  }
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
        t.par.clfixed_x, t.par.coalth, {}, t.par.verb,
        t.par.recolor_unionfind, t.par.recolor_reduce, t.par.recolor_grid, m);
  }
  if (csv_out && sem.Nested()) {
    auto plic = t.solver->GetPlic();
    CalcTraj<M>(
        m, plic.layers, plic.vfcu, &t.fccl, plic.vfcim, t.fc_zero,
        t.fc_zerov, t.column_names, t.row_colors, t.table);
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
  const bool isroot = true;

  ArgumentParser parser("Sharpens the image using PLIC advection", isroot);
  parser.AddSwitch({"--verbose", "-v"}).Help("Report steps");
  parser.AddVariable<double>("--cfl", 0.1)
      .Help("CFL number for advection, valid values between 0 and 1");
  parser.AddVariable<std::string>("--vtk_out_march")
      .Help("Path to output VTK with surface from marching cubes");
  parser.AddVariable<std::string>("--vtk_out").Help(
      "Path to output VTK with piecewise linear surface");
  parser.AddVariable<std::string>("--csv_out").Help(
      "Path to output CSV with centroids of connected components");
  parser.AddVariable<int>("--steps", 5).Help("Number of sharpening steps");

  parser.AddVariable<std::string>("hdf_in").Help(
      "Path to input image as HDF5 array of floats between 0 and 1");
  parser.AddVariable<std::string>("hdf_out").Help(
      "Path to output image, can be the same as HDF_IN");

  auto args = parser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  const auto shape = Hdf<M>::GetShape(args.String["hdf_in"]);
  const int nx = shape[2];
  const int ny = shape[1];
  const int nz = shape[0];

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
    while (nx % bs == 0 && ny % bs == 0 && bs <= 32) {
      bs *= 2;
    }
    bs /= 2;
    bsz = bs;
  }
  int bsy = bs;

  std::stringstream conf;
  conf << R"EOF(
# ranks
set int px 1
set int py 1
set int pz 1

set int hl 2

set int verbose_time 0
set int verbose_stages 0
set string backend cubismnc
set int loc_maxcomm 16
)EOF";

  conf << "set int bx " << nx / bs << '\n';
  conf << "set int by " << ny / bsy << '\n';
  conf << "set int bz " << nz / bsz << '\n';

  conf << "set int bsx " << bs << '\n';
  conf << "set int bsy " << bsy << '\n';
  conf << "set int bsz " << bsz << '\n';

  conf << "set int dim " << (nz == 1 ? 2 : 3) << '\n';

  conf << "set int steps " << args.Int["steps"] << '\n';

  conf << "set string hdf_in " << args.String["hdf_in"] << '\n';
  conf << "set string hdf_out " << args.String["hdf_out"] << '\n';
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

  MpiWrapper mpi(&argc, &argv);
  return RunMpiBasicString<M>(mpi, Run, conf.str());
}
