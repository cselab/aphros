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

// fcu: quantity to advect
// mec: boundary conditions
// vel: advection velocity
// dt: time step
void Advection0(
    FieldCell<Scal>& fcu, const MapCondFace& mec, const FieldEmbed<Scal>& fev,
    Scal dt, const Embed<M>& eb) {
  const auto& m = eb.GetMesh();
  const auto feu = eb.InterpolateUpwind(fcu, fev.GetFieldFace(), mec, 0, 1.);
  // Compute flux.
  FieldEmbed<Scal> fevu(m, 0);
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
      sum += fevu[m.GetFace(c, q)] * m.GetOutwardFactor(c, q);
    }
    fcu[c] -= sum * dt / eb.GetVolume(c);
  }
}

void Run(M& m, Vars& var) {
  auto sem = m.GetSem();
  struct {
    FieldNode<Scal> fnl;
    FieldCell<Scal> fcu;
    FieldEmbed<Scal> fev;
    MapCondFace mfc; // boundary conditions
    Scal t = 0;
    Vect vel;
    Scal dt;
    std::unique_ptr<Embed<M>> eb_;
  } * ctx(sem);
  auto& fnl = ctx->fnl;
  auto& fcu = ctx->fcu;
  auto& fev = ctx->fev;
  auto& mfc = ctx->mfc;
  auto& t = ctx->t;
  auto& dt = ctx->dt;
  auto& vel = ctx->vel;

  if (sem("init")) {
    vel = Vect(var.Vect["vel"]);
    dt = var.Double["cfl"] * m.GetCellSize()[0] / vel.norm1();
    fcu.Reinit(m, 0);
    for (auto c : m.AllCells()) {
      auto x = m.GetCenter(c);
      fcu[c] = ((Vect(0.4, 0.4, 0) - x).norminf() < 0.2);
    }
    ctx->eb_.reset(new Embed<M>(m));
    fnl.Reinit(m, 0);
    for (auto n : m.AllNodes()) {
      auto x = m.GetNode(n);
      fnl[n] = (0.4 - (Vect(0.5, 0.5, 0) - x).norm());
    }
    m.Dump(&fcu, "u");
  }
  if (sem.Nested("eb-init")) {
    ctx->eb_->Init(fnl);
  }
  if (sem("flux")) {
    auto& eb = *ctx->eb_;
    fev.Reinit(m, 0);
    for (auto f : eb.Faces()) {
      fev[f] = eb.GetSurface(f).dot(vel);
    }
    for (auto c : eb.CFaces()) {
      fev[c] = eb.GetSurface(c).dot(vel);
    }
  }
  if (sem.Nested("eb-dump")) {
    ctx->eb_->DumpPoly();
  }
  sem.LoopBegin();
  if (sem("step")) {
    switch (var.Int["example"]) {
      case 0:
        Advection0(fcu, mfc, fev, dt, *ctx->eb_);
        break;
      default:
        throw std::runtime_error(
            "Unknown example=" + std::to_string(var.Int["example"]));
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

set double tmax 0.2
set vect vel 1 1 0
set double cfl 0.5

set int example 0
)EOF";

  if (argc > 1) {
    const int example = atoi(argv[1]);
    conf += "\nset int example " + std::to_string(example);
  }

  return RunMpiBasic<M>(argc, argv, Run, conf);
}
