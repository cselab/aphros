#undef NDEBUG
#include <cassert>
#include <iostream>
#include <string>
#include <memory>

#include "distr/distrbasic.h"
#include "solver/embed.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;

// fcu: quantity to advect
// mec: boundary conditions
// vel: advection velocity
// dt: time step
void Advection(FieldCell<Scal>& fcu, const MapCondFace& mec,
               Vect vel, Scal dt, const Embed<M>& eb) {
  const auto& m = eb.GetMesh();
  const auto feu = eb.Interpolate(fcu, mec, 0, 0.);
  // Compute flux through all faces, zero in regular cells.
  FieldEmbed<Scal> fevu(m, 0);
  for (auto f : eb.Faces()) {
    fevu[f] = feu[f] * vel.dot(eb.GetSurface(f));
  }
  for (auto c : eb.CFaces()) {
    fevu[c] = feu[c] * vel.dot(eb.GetSurface(c));
  }
  // Advance in time.
  for (auto c : eb.Cells()) {
    Scal sum = fevu[c];
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
    MapCondFace mfc; // face conditions
    Scal t = 0;
    std::unique_ptr<Embed<M>> eb_;
  } * ctx(sem);

  auto& fcu = ctx->fcu;
  const Scal tmax = 0.3;
  const Vect vel(0.5, 0.3, 0.1);
  const Scal cfl = 0.5;

  if (sem("init")) {
    fcu.Reinit(m, 0);
    for (auto c : m.Cells()) {
      fcu[c] = (m.GetCenter(c).dist(Vect(0.5, 0.5, 0.5)) < 0.2);
    }
    ctx->eb_.reset(new Embed<M>(m));
    ctx->fnl = InitEmbed(m, var, m.IsRoot());
  }
  if (sem.Nested("eb-init")) {
    ctx->eb_->Init(ctx->fnl);
  }
  if (sem.Nested("eb-dump")) {
    ctx->eb_->DumpPoly();
  }
  sem.LoopBegin();
  if (sem("step")) {
    const Scal dt = cfl * m.GetCellSize()[0] / vel.norm1();
    Advection(fcu, ctx->mfc, vel, dt, *ctx->eb_);
    ctx->t += dt;
  }
  if (sem("checkloop")) {
    if (ctx->t >= tmax) {
      sem.LoopBreak();
    }
  }
  if (sem("dump")) {
    m.Dump(&fcu, "u");
  }
  if (sem()) {
  }
  sem.LoopEnd();
}

int main(int argc, const char** argv) {
  std::string conf = R"EOF(
set int bx 1
set int by 1
set int bz 1

set int bsx 32
set int bsy 32
set int bsz 32

set int px 1
set int py 1
set int pz 1

set string eb_init sphere
set vect eb_sphere_c 0.5 0.5 0.5
set vect eb_sphere_r 0.4 0.4 0.4
set double eb_sphere_angle 0
set int eb_init_inverse 1
)EOF";

  return RunMpiBasic<M>(argc, argv, Run, conf);
}
