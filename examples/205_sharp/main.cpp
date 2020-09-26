// Created by Petr Karnakov on 25.09.2020
// Copyright 2020 ETH Zurich

#undef NDEBUG
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
  std::string vtk_out = var.String("vtk_out", "");
  if (vtk_out.length() && sem.Nested()) {
    t.solver->DumpInterface(vtk_out);
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

  std::string hdf_in;
  int steps;
  int nx, ny, nz;
  std::string hdf_out;

  const std::set<std::string> novalue_args = {
      "-v",
      "--verbose",
      "-h",
      "--help",
  };
  const auto args = ParseArgs(argc, argv, novalue_args);
  const std::map<std::string, std::string> oargs = args.first; // optional
  const std::vector<std::string> pargs = args.second; // positional

  auto print_usage = [&argv, isroot](bool full) {
    auto& s = std::cerr;
    if (isroot) {
      s << "usage: " << argv[0]
        << " [-h|--help]"
           " [-v|--verbose]"
           " [--vtk_out VTK_OUT]"
           " [--cfl CFL]"
           " HDF_IN NX NY NZ STEPS HDF_OUT"
           "\n";
      if (full) {
        s << 
          "Sharpens the image using PLIC advection.\n"
          "HDF_IN: path to input image as HDF5 array of floats between 0 and 1 and shape (1,NZ,NY,NX)\n"
          "HDF_OUT: path to output image\n"
          "VTK_OUT: path to output vtk with piecewise linear surface\n"
          "STEPS: number of sharpening steps\n"
          "CFL: CFL number for advection, valid values between 0 and 1\n"
          ;
      }
    }
  };

  if (oargs.count("-h") || oargs.count("--help")) {
    print_usage(true);
    return 0;
  }

  if (pargs.size() != 6) {
    if (isroot) {
      std::cerr << "invalid number of arguments: " << pargs.size() << "\n";
    }
    print_usage(false);
    return 1;
  }

  auto known_args = novalue_args;
  known_args.insert("--vtk_out");
  known_args.insert("--cfl");
  auto check_known_args = [&oargs, isroot, known_args]() {
    bool pass = true;
    for (auto p : oargs) {
      if (!known_args.count(p.first)) {
        pass = false;
        if (isroot) {
          std::cerr << "unrecognized option: " << p.first << '\n';
        }
      }
    }
    return pass;
  };
  if (!check_known_args()) {
    print_usage(false);
    return 1;
  }

  int i = 0;
  hdf_in = pargs[i++];
  nx = GetInt(pargs[i++]);
  ny = GetInt(pargs[i++]);
  nz = GetInt(pargs[i++]);
  steps = GetInt(pargs[i++]);
  hdf_out = pargs[i++];

  /*
  auto checkdiv = [](int n, int b, std::string name) {
    if (n % b) {
      std::stringstream s;
      s << name << "=" << n << " not divisible by " << b;
      throw std::runtime_error(s.str());
    }
  };

  int bs = 8;
  int bsy = bs;
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
  */

  // XXX: adhoc for backend=local
  int bs = nx;
  int bsy = ny;
  int bsz = nz;

  const bool verbose_steps = oargs.count("-v") || oargs.count("--verbose");

  std::stringstream conf;
  conf << R"EOF(
# ranks
set int px 1
set int py 1
set int pz 1

set int hl 2

set int verbose_time 0
set int verbose_stages 0
set string backend local
set int loc_maxcomm 16
)EOF";

  conf << "set int bx " << nx / bs << '\n';
  conf << "set int by " << ny / bsy << '\n';
  conf << "set int bz " << nz / bsz << '\n';

  conf << "set int bsx " << bs << '\n';
  conf << "set int bsy " << bsy << '\n';
  conf << "set int bsz " << bsz << '\n';

  conf << "set int dim " << (nz == 1 ? 2 : 3) << '\n';

  conf << "set int steps " << steps << '\n';

  conf << "set string hdf_in " << hdf_in << '\n';
  conf << "set string hdf_out " << hdf_out << '\n';
  if (oargs.count("--vtk_out")) {
    conf << "set string vtk_out " << oargs.at("--vtk_out") << '\n';
  }
  {
    double cfl = 0.5;
    if (oargs.count("--cfl")) {
      cfl = GetDouble(oargs.at("--cfl"));
    }
    conf << "set double cfl " << cfl << '\n';
  }
  conf << "set int VERBOSE " << (verbose_steps ? 1 : 0) << '\n';

  return RunMpiBasic<M>(argc, argv, Run, conf.str());
}
