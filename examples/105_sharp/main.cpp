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

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;

void Run(M& m, Vars& var) {
  const bool verbose = var.Int["VERBOSE"];
  using Scal = typename M::Scal;
  auto sem = m.GetSem();
  struct {
    FieldNode<Scal> fnl;
    FieldCell<Scal> fcvx;
    std::unique_ptr<Vof<M>> solver;

    FieldCell<Scal> fcu;
    FieldEmbed<Scal> fe_flux;
    FieldCell<Scal> fc_src;

    // boundary conditions for advection (empty)
    MapEmbed<BCondAdvection<Scal>> bc;
  } * ctx(sem);
  auto& t = *ctx;
  if (sem.Nested()) {
    Hdf<M>::Read(t.fcu, var.String["hdf_in"], m);
  }
  if (sem("ctor")) {
    t.fe_flux.Reinit(m, 0);
    t.fc_src.Reinit(m, 0);

    Vof<M>::Par par;
    par.clipth = 1e-10;
    par.sharpen = true;
    par.sharpen_cfl = var.Double["cfl"];
    par.dim = var.Int["dim"];
    const FieldCell<Scal> fccl(m, 0);
    t.solver.reset(new Vof<M>(
        m, m, t.fcu, fccl, t.bc, &t.fe_flux, &t.fc_src, 0., 1., par));
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
  if (sem()) { // XXX empty stage
  }
}

int GetInt(std::string s) {
  int r;
  std::stringstream t(s);
  t >> r;
  return r;
}

double GetDouble(std::string s) {
  double r;
  std::stringstream t(s);
  t >> r;
  return r;
}

int main(int argc, const char** argv) {
  const bool isroot = true;

  ArgumentParser parser("Sharpens the image using PLIC advection", isroot);
  parser.AddSwitch({"--verbose", "-v"}).Help("Report steps");
  parser.AddVariable<double>("--cfl", 0.5)
      .Help("CFL number for advection, valid values between 0 and 1");
  parser.AddVariable<std::string>("--vtk_out_march")
      .Help("Path to output VTK with surface from marching cubes");
  parser.AddVariable<std::string>("--vtk_out").Help(
      "Path to output VTK with piecewise linear surface");

  parser.AddVariable<std::string>("hdf_in").Help(
      "Path to input image as HDF5 array of floats between 0 and 1");
  parser.AddVariable<int>("steps", 1).Help("Number of sharpening steps");
  parser.AddVariable<std::string>("hdf_out").Help("Path to output image");

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
  conf << "set double cfl " << args.Double["cfl"] << '\n';
  conf << "set int VERBOSE " << args.Int["verbose"] << '\n';

  return RunMpiBasic<M>(argc, argv, Run, conf.str());
}
