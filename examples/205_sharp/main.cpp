// Created by Petr Karnakov on 05.07.2020
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
    par.clipth = 0;
    par.sharpen = true;
    par.sharpen_cfl = 0.1;
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
  if (sem()) { // XXX empty stage
  }
}

int GetInt(std::string s) {
  int r;
  std::stringstream t(s);
  t >> r;
  return r;
}

int main(int argc, const char** argv) {
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

  if (pargs.size() != 6 || oargs.count("-h") || oargs.count("--help")) {
    std::cerr << "usage: " << argv[0]
              << " [-h|--help]"
                 " [-v|--verbose]"
                 " HDF_IN NX NY NZ STEPS HDF_OUT"
              << std::endl;
    return 1;
  } else {
    int i = 0;
    hdf_in = pargs[i++];
    nx = GetInt(pargs[i++]);
    ny = GetInt(pargs[i++]);
    nz = GetInt(pargs[i++]);
    steps = GetInt(pargs[i++]);
    hdf_out = pargs[i++];
  }

  auto checkdiv = [](int n, int b, std::string name) {
    if (n % b) {
      std::cerr << name << "=" << n << " not divisible by " << b << '\n';
      std::terminate();
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
)EOF";

  conf << "set int bx " << nx / bs << '\n';
  conf << "set int by " << ny / bs << '\n';
  conf << "set int bz " << nz / bsz << '\n';

  conf << "set int bsx " << bs << '\n';
  conf << "set int bsy " << bs << '\n';
  conf << "set int bsz " << bsz << '\n';

  conf << "set int steps " << steps << '\n';

  conf << "set string hdf_in " << hdf_in << '\n';
  conf << "set string hdf_out " << hdf_out << '\n';
  conf << "set int VERBOSE " << (verbose_steps ? 1 : 0) << '\n';

  return RunMpiBasic<M>(argc, argv, Run, conf.str());
}
