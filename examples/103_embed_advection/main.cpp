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
    FieldCell<Scal>& fcu, const MapCondFace& mec, Vect vel, Scal dt,
    const Embed<M>& eb) {
  const auto& m = eb.GetMesh();
  const auto feu = eb.Interpolate(fcu, mec, 0, 0.);
  // Compute flux.
  FieldEmbed<Scal> fevu(m, 0); // quantity flux
  for (auto f : eb.Faces()) {
    fevu[f] = feu[f] * vel.dot(eb.GetSurface(f));
  }
  for (auto c : eb.CFaces()) {
    fevu[c] = feu[c] * vel.dot(eb.GetSurface(c));
  }
  // Advance in time.
  for (auto c : eb.Cells()) {
    Scal sum = fevu[c]; // sum of fluxes
    for (auto q : eb.Nci(c)) {
      sum += fevu[m.GetFace(c, q)] * m.GetOutwardFactor(c, q);
    }
    fcu[c] += sum * dt / eb.GetVolume(c);
  }
}

void Run(M& m, Vars& var) {
  auto sem = m.GetSem();
  struct {
    FieldNode<Scal> fnl;
    FieldCell<Scal> fcu;
    MapCondFace mfc; // boundary conditions
    Scal t = 0;
    std::unique_ptr<Embed<M>> eb_;
  } * ctx(sem);

  if (sem("init")) {
    ctx->fcu.Reinit(m, 0);
    for (auto c : m.AllCells()) {
      auto x = m.GetCenter(c);
      ctx->fcu[c] = ((Vect(0.5, 0.5, 0) - x).norminf() < 0.2);
    }
    ctx->eb_.reset(new Embed<M>(m));
    ctx->fnl.Reinit(m, 0);
    for (auto n : m.AllNodes()) {
      auto x = m.GetNode(n);
      ctx->fnl[n] = (0.4 - (Vect(0.5, 0.5, 0) - x).norm());
    }
  }
  if (sem.Nested("eb-init")) {
    ctx->eb_->Init(ctx->fnl);
  }
  if (sem.Nested("eb-dump")) {
    ctx->eb_->DumpPoly();
  }
  sem.LoopBegin();
  if (sem("step")) {
    const Vect vel(var.Vect["vel"]);
    const Scal dt = var.Double["cfl"] * m.GetCellSize()[0] / vel.norm1();
    switch (var.Int["example"]) {
      case 0:
        Advection0(ctx->fcu, ctx->mfc, vel, dt, *ctx->eb_);
        break;
      default:
        throw std::runtime_error(
            "Unknown example=" + std::to_string(var.Int["example"]));
    }
    ctx->t += dt;
  }
  if (sem("loopbreak")) {
    if (ctx->t >= var.Double["tmax"]) {
      sem.LoopBreak();
    }
  }
  if (sem("dump")) {
    m.Dump(&ctx->fcu, "u");
  }
  sem.LoopEnd();
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

set double tmax 1
set vect vel 0.5 0.3 0
set double cfl 0.5

set int example 0
)EOF";

  if (argc > 1) {
    const int example = atoi(argv[1]);
    conf += "\nset int example " + std::to_string(example);
  }

  return RunMpiBasic<M>(argc, argv, Run, conf);
}
