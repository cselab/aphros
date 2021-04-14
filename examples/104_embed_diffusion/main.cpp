// Created by Petr Karnakov on 29.01.2020
// Copyright 2020 ETH Zurich

#include <cassert>
#include <cstdio>
#include <iostream>
#include <memory>
#include <string>

#include "distr/distrbasic.h"
#include "parse/argparse.h"
#include "solver/approx_eb.h"
#include "solver/embed.h"

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;

// Diffusion solver.
// fcu: quantity to advect
// mebc: boundary conditions
// diff: diffusion coefficient
// dt: time step
void Diffusion0(
    FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc, Scal diff, Scal dt,
    const Embed<M>& eb) {
  const auto feg = UEmbed<M>::Gradient(fcu, mebc, eb);
  // Compute flux.
  FieldEmbed<Scal> fed(eb, 0);
  for (auto f : eb.Faces()) {
    fed[f] = feg[f] * diff * eb.GetArea(f);
  }
  for (auto c : eb.CFaces()) {
    fed[c] = feg[c] * diff * eb.GetArea(c);
  }
  // Advance in time.
  for (auto c : eb.Cells()) {
    Scal sum = fed[c]; // sum of fluxes
    for (auto q : eb.Nci(c)) {
      sum += fed[eb.GetFace(c, q)] * eb.GetOutwardFactor(c, q);
    }
    fcu[c] += sum * dt / eb.GetVolume(c);
  }
}

// Diffusion solver with redistribution to neighbor cells.
void Diffusion1(
    FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc, Scal diff, Scal dt,
    const Embed<M>& eb) {
  const auto feg = UEmbed<M>::Gradient(fcu, mebc, eb);
  // Compute flux.
  FieldEmbed<Scal> fed(eb, 0);
  for (auto f : eb.Faces()) {
    fed[f] = feg[f] * diff * eb.GetArea(f);
  }
  for (auto c : eb.CFaces()) {
    fed[c] = feg[c] * diff * eb.GetArea(c);
  }
  // Compute the change at one time step.
  FieldCell<Scal> fct(eb, 0);
  for (auto c : eb.Cells()) {
    Scal sum = fed[c];
    for (auto q : eb.Nci(c)) {
      sum += fed[eb.GetFace(c, q)] * eb.GetOutwardFactor(c, q);
    }
    fct[c] = sum * dt;
  }
  fct = UEmbed<M>::RedistributeCutCells(fct, eb);
  // Advance in time.
  for (auto c : eb.Cells()) {
    fcu[c] += fct[c] / eb.GetVolume(c);
  }
}

// Diffusion solver with redistribution to neighbor cells
// and LoopFaces() to avoid code duplicatoin.
void Diffusion2(
    FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc, Scal diff, Scal dt,
    const Embed<M>& eb) {
  const auto feg = UEmbed<M>::Gradient(fcu, mebc, eb);
  // Compute flux.
  FieldEmbed<Scal> fed(eb, 0);
  eb.LoopFaces([&](auto cf) { // lambda-function applied to faces and
                              // embedded faces
    fed[cf] = feg[cf] * diff * eb.GetArea(cf);
  });
  // Compute the change at one time step.
  FieldCell<Scal> fct(eb, 0);
  for (auto c : eb.Cells()) {
    Scal sum = fed[c];
    for (auto q : eb.Nci(c)) {
      sum += fed[eb.GetFace(c, q)] * eb.GetOutwardFactor(c, q);
    }
    fct[c] = sum * dt;
  }
  fct = UEmbed<M>::RedistributeCutCells(fct, eb);
  // Advance in time.
  for (auto c : eb.Cells()) {
    fcu[c] += fct[c] / eb.GetVolume(c);
  }
}

void Run(M& m, Vars& var) {
  auto sem = m.GetSem();
  struct {
    FieldNode<Scal> fnl;
    FieldCell<Scal> fcu;
    MapEmbed<BCond<Scal>> mebc;
    Scal t = 0;
    Scal diff;
    Scal dt;
    std::unique_ptr<Embed<M>> eb_;
  } * ctx(sem);
  auto& fnl = ctx->fnl;
  auto& fcu = ctx->fcu;
  auto& mebc = ctx->mebc;
  auto& t = ctx->t;
  auto& dt = ctx->dt;
  auto& diff = ctx->diff;

  if (sem("init")) {
    diff = var.Double["diff"];
    dt = var.Double["cfl"] * sqr(m.GetCellSize()[0]) * diff;
    // initial field
    fcu.Reinit(m, 0);
    for (auto c : m.AllCells()) {
      auto x = m.GetCenter(c);
      fcu[c] = ((Vect(0.5, 0.5, 0) - x).norminf() < 0.1);
    }
    // level-set for embedded boundaries
    fnl.Reinit(m, 0);
    for (auto n : m.AllNodes()) {
      auto x = m.GetNode(n);
      fnl[n] = (0.4 - (Vect(0.5, 0.5, 0) - Vect(x[0], x[1], 0.)).norm());
    }
    ctx->eb_.reset(new Embed<M>(m, var.Double["gradlim"]));
    m.Dump(&fcu, "u");
  }
  if (sem.Nested("eb-init")) {
    ctx->eb_->Init(fnl);
  }
  if (sem.Nested("eb-dump")) {
    ctx->eb_->DumpPlaneSection(Vect(0., 0., 1e-3), Vect(0., 0., 1.));
  }
  if (sem("init2")) {
    auto& eb = *ctx->eb_;
    // boundary conditions
    for (auto c : eb.SuCFaces()) {
      auto& bc = mebc[c];
      bc.type = BCondType::dirichlet;
      bc.val = 1;
    }
  }
  sem.LoopBegin();
  if (sem("step")) {
    if (m.IsRoot()) {
      std::cout << "t=" << t << std::endl;
    }
    switch (var.Int["case"]) {
      case 0:
        Diffusion0(fcu, mebc, diff, dt, *ctx->eb_);
        break;
      case 1:
        Diffusion1(fcu, mebc, diff, dt, *ctx->eb_);
        break;
      case 2:
        Diffusion2(fcu, mebc, diff, dt, *ctx->eb_);
        break;
      default:
        throw std::runtime_error(
            "Unknown case=" + std::to_string(var.Int["case"]));
    }
    t += dt;
    m.Comm(&fcu);
  }
  if (sem("loopbreak")) {
    if (t >= var.Double["tmax"]) {
      sem.LoopBreak();
      m.Dump(&fcu, "u");
    }
  }
  sem.LoopEnd();
  if (sem()) {
  }
}

int main(int argc, const char** argv) {
  std::string conf = R"EOF(
# domain size in blocks
set int bx 1
set int by 1
set int bz 1

# block size in cells
set int bsx 32
set int bsy 32
set int bsz 1

set string dumpformat plain

set double tmax 0.002
set double diff 1
set double cfl 0.25

)EOF";

  ArgumentParser parser("Diffusion example with embedded boundaries");
  parser.AddVariable<int>("case", 0).Help("Case to run, choices: 0, 1, 2");
  parser.AddVariable<double>("gradlim", 0).Help("Value of gradient limiter");
  auto args = parser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  conf += "\nset int case " + args.Int.GetStr("case");
  conf += "\nset double gradlim " + args.Double.GetStr("gradlim");

  MpiWrapper mpi(&argc, &argv);
  return RunMpiBasicString<M>(mpi, Run, conf);
}
