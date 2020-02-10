// Created by Petr Karnakov on 29.01.2020
// Copyright 2020 ETH Zurich

#include <cassert>
#include <cstdio>
#include <iostream>
#include <memory>
#include <string>

#include "distr/distrbasic.h"
#include "solver/embed.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;

// Diffusion solver.
// fcu: quantity to advect
// mec: boundary conditions
// diff: diffusion coefficient
// dt: time step
void Diffusion0(
    FieldCell<Scal>& fcu, const MapCondFace& mec, Scal diff, Scal dt,
    const Embed<M>& eb) {
  const auto& m = eb.GetMesh();
  const auto feg = eb.Gradient(fcu, mec, 0, 1.);
  // Compute flux.
  FieldEmbed<Scal> fed(m, 0);
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
      sum += fed[m.GetFace(c, q)] * m.GetOutwardFactor(c, q);
    }
    fcu[c] += sum * dt / eb.GetVolume(c);
  }
}

void Diffusion1(
    FieldCell<Scal>& fcu, const MapCondFace& mec, Scal diff, Scal dt,
    const Embed<M>& eb) {
  const auto& m = eb.GetMesh();
  const auto feg = eb.Gradient(fcu, mec, 0, 1.);
  // Compute flux.
  FieldEmbed<Scal> fed(m, 0);
  for (auto f : eb.Faces()) {
    fed[f] = feg[f] * diff * eb.GetArea(f);
  }
  for (auto c : eb.CFaces()) {
    fed[c] = feg[c] * diff * eb.GetArea(c);
  }
  // Compute the change at one time step.
  FieldCell<Scal> fct(m, 0);
  for (auto c : eb.Cells()) {
    Scal sum = fed[c];
    for (auto q : eb.Nci(c)) {
      sum += fed[m.GetFace(c, q)] * m.GetOutwardFactor(c, q);
    }
    fct[c] = sum * dt;
  }
  fct = eb.RedistributeCutCells(fct);
  // Advance in time.
  for (auto c : eb.Cells()) {
    fcu[c] += fct[c] / eb.GetVolume(c);
  }
}

void Diffusion2(
    FieldCell<Scal>& fcu, const MapCondFace& mec, Scal diff, Scal dt,
    const Embed<M>& eb) {
  const auto& m = eb.GetMesh();
  const auto feg = eb.GradientBilinear(fcu, mec, 0, 1.);
  // Compute flux.
  FieldEmbed<Scal> fed(m, 0);
  for (auto f : eb.Faces()) {
    fed[f] = feg[f] * diff * eb.GetArea(f);
  }
  for (auto c : eb.CFaces()) {
    fed[c] = feg[c] * diff * eb.GetArea(c);
  }
  // Compute the change at one time step.
  FieldCell<Scal> fct(m, 0);
  for (auto c : eb.Cells()) {
    Scal sum = fed[c];
    for (auto q : eb.Nci(c)) {
      sum += fed[m.GetFace(c, q)] * m.GetOutwardFactor(c, q);
    }
    fct[c] = sum * dt;
  }
  fct = eb.RedistributeCutCells(fct);
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
    MapCondFace mfc; // boundary conditions
    Scal t = 0;
    Scal diff;
    Scal dt;
    std::unique_ptr<Embed<M>> eb_;
  } * ctx(sem);
  auto& fnl = ctx->fnl;
  auto& fcu = ctx->fcu;
  auto& mfc = ctx->mfc;
  auto& t = ctx->t;
  auto& dt = ctx->dt;
  auto& diff = ctx->diff;

  if (sem("init")) {
    diff = var.Double["diff"];
    dt = var.Double["cfl"] * sqr(m.GetCellSize()[0]) * diff;
    fcu.Reinit(m, 0);
    for (auto c : m.AllCells()) {
      auto x = m.GetCenter(c);
      fcu[c] = ((Vect(0.5, 0.5, 0) - x).norminf() < 0.1);
    }
    ctx->eb_.reset(new Embed<M>(m, var.Double["gradlim"]));
    fnl.Reinit(m, 0);
    for (auto n : m.AllNodes()) {
      auto x = m.GetNode(n);
      fnl[n] = (0.4 - (Vect(0.5, 0.5, 0) - Vect(x[0], x[1], 0.)).norm());
    }
    m.Dump(&fcu, "u");
  }
  if (sem.Nested("eb-init")) {
    ctx->eb_->Init(fnl);
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
        Diffusion0(fcu, mfc, diff, dt, *ctx->eb_);
        break;
      case 1:
        Diffusion1(fcu, mfc, diff, dt, *ctx->eb_);
        break;
      case 2:
        Diffusion2(fcu, mfc, diff, dt, *ctx->eb_);
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
set double cfl 0.2

set double gradlim 0

set int case 0
)EOF";

  if (argc > 1) {
    const int case_ = atoi(argv[1]);
    conf += "\nset int case " + std::to_string(case_);
  }

  if (argc > 2) {
    const double gradlim = atof(argv[2]);
    conf += "\nset double gradlim " + std::to_string(gradlim);
  }

  return RunMpiBasic<M>(argc, argv, Run, conf);
}
