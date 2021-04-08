// Created by Petr Karnakov on 29.01.2020
// Copyright 2020 ETH Zurich

#include <cassert>
#include <cstdio>
#include <iostream>
#include <memory>
#include <string>

#include <distr/distrbasic.h>
#include <parse/argparse.h>
#include <solver/approx_eb.h>
#include <solver/embed.h>

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;

// Advection solver.
// fcu: quantity to advect
// mebc: boundary conditions
// vel: advection velocity
// dt: time step
void Advection0(
    FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc,
    const FieldEmbed<Scal>& fev, Scal dt, const Embed<M>& eb) {
  const auto fcg =
      UEmbed<M>::GradientLinearFit(UEmbed<M>::Interpolate(fcu, mebc, eb), eb);
  const auto feu =
      UEmbed<M>::InterpolateUpwind(fcu, mebc, ConvSc::fou, fcg, fev, eb);
  // Compute flux.
  FieldEmbed<Scal> fevu(eb, 0);
  for (auto f : eb.Faces()) {
    fevu[f] = feu[f] * fev[f];
  }
  for (auto c : eb.CFaces()) {
    fevu[c] = feu[c] * fev[c];
  }
  // Advance in time.
  for (auto c : eb.Cells()) {
    Scal sum = fevu[c]; // sum of fluxes
    for (auto q : eb.Nci(c)) {
      sum += fevu[eb.GetFace(c, q)] * eb.GetOutwardFactor(c, q);
    }
    fcu[c] -= sum * dt / eb.GetVolume(c);
  }
}

// Advection solver with redistribution from cut cells
// and first order upwind scheme.
void Advection1(
    FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc,
    const FieldEmbed<Scal>& fev, Scal dt, const Embed<M>& eb) {
  const auto fcg =
      UEmbed<M>::GradientLinearFit(UEmbed<M>::Interpolate(fcu, mebc, eb), eb);
  const auto feu =
      UEmbed<M>::InterpolateUpwind(fcu, mebc, ConvSc::fou, fcg, fev, eb);
  // Compute flux.
  FieldEmbed<Scal> fevu(eb, 0);
  for (auto f : eb.Faces()) {
    fevu[f] = feu[f] * fev[f];
  }
  for (auto c : eb.CFaces()) {
    fevu[c] = feu[c] * fev[c];
  }
  // Compute the change at one time step.
  FieldCell<Scal> fct(eb, 0);
  for (auto c : eb.Cells()) {
    Scal sum = fevu[c];
    for (auto q : eb.Nci(c)) {
      sum += fevu[eb.GetFace(c, q)] * eb.GetOutwardFactor(c, q);
    }
    fct[c] = -sum * dt;
  }
  fct = UEmbed<M>::RedistributeCutCells(fct, eb);
  // Advance in time.
  for (auto c : eb.Cells()) {
    fcu[c] += fct[c] / eb.GetVolume(c);
  }
}

// Advection solver with redistribution from cut cells
// and second order upwind scheme.
void Advection2(
    FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc,
    const FieldEmbed<Scal>& fev, Scal dt, const Embed<M>& eb) {
  const auto fcg =
      UEmbed<M>::GradientLinearFit(UEmbed<M>::Interpolate(fcu, mebc, eb), eb);
  auto feu = UEmbed<M>::InterpolateUpwind(fcu, mebc, ConvSc::sou, fcg, fev, eb);
  // Compute flux.
  FieldEmbed<Scal> fevu(eb, 0);
  for (auto f : eb.Faces()) {
    fevu[f] = feu[f] * fev[f];
  }
  for (auto c : eb.CFaces()) {
    fevu[c] = feu[c] * fev[c];
  }
  // Compute the change at one time step.
  FieldCell<Scal> fct(eb, 0);
  for (auto c : eb.Cells()) {
    Scal sum = fevu[c];
    for (auto q : eb.Nci(c)) {
      sum += fevu[eb.GetFace(c, q)] * eb.GetOutwardFactor(c, q);
    }
    fct[c] = -sum * dt;
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
    FieldEmbed<Scal> fev;
    MapEmbed<BCond<Scal>> mebc;
    Scal t = 0;
    Vect vel;
    Scal dt;
    std::unique_ptr<Embed<M>> eb_;
  } * ctx(sem);
  auto& fnl = ctx->fnl;
  auto& fcu = ctx->fcu;
  auto& fev = ctx->fev;
  auto& mebc = ctx->mebc;
  auto& t = ctx->t;
  auto& dt = ctx->dt;
  auto& vel = ctx->vel;

  if (sem("init")) {
    vel = Vect(var.Vect["vel"]);
    dt = var.Double["cfl"] * m.GetCellSize()[0] / vel.norm1();
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
    ctx->eb_.reset(new Embed<M>(m));
    m.Dump(&fcu, "u");
  }
  if (sem.Nested("eb-init")) {
    ctx->eb_->Init(fnl);
  }
  if (sem("init2")) {
    auto& eb = *ctx->eb_;
    // boundary conditions
    for (auto c : eb.SuCFaces()) {
      auto& bc = mebc[c];
      bc.type = BCondType::dirichlet;
      bc.val = 1;
    }
    // flux
    fev.Reinit(m, 0);
    for (auto f : eb.Faces()) {
      fev[f] = eb.GetSurface(f).dot(vel);
    }
    for (auto c : eb.CFaces()) {
      fev[c] = eb.GetSurface(c).dot(vel);
    }
  }
  if (sem.Nested("eb-dump")) {
    ctx->eb_->DumpPlaneSection(Vect(0., 0., 1e-3), Vect(0., 0., 1.));
  }
  sem.LoopBegin();
  if (sem("step")) {
    if (m.IsRoot()) {
      std::cout << "t=" << t << std::endl;
    }
    switch (var.Int["case"]) {
      case 0:
        Advection0(fcu, mebc, fev, dt, *ctx->eb_);
        break;
      case 1:
        Advection1(fcu, mebc, fev, dt, *ctx->eb_);
        break;
      case 2:
        Advection2(fcu, mebc, fev, dt, *ctx->eb_);
        break;
      default:
        throw std::runtime_error(
            "Unknown example=" + std::to_string(var.Int["case"]));
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

set double tmax 0.1
set vect vel 2 1 0
set double cfl 0.25
)EOF";

  ArgumentParser parser("Advection example with embedded boundaries");
  parser.AddVariable<int>("case", 0).Help("Case to run, choices: 0, 1, 2");
  auto args = parser.ParseArgs(argc, argv);
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  conf += "\nset int case " + args.Int.GetStr("case");

  MpiWrapper mpi(&argc, &argv);
  return RunMpiBasicString<M>(mpi, Run, conf);
}
