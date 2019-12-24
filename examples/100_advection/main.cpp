#undef NDEBUG
#include <cassert>
#include <iostream>
#include <string>
#include <memory>

#include "distr/distrbasic.h"
#include "solver/vof.h"

using namespace solver;

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using namespace solver;

struct State {};

void Run(M& m, State&, Vars&) {
  auto sem = m.GetSem();

  struct {
    std::unique_ptr<Vof<M>> as; // advection solver
    FieldCell<Scal> fc_src; // volume source
    FieldFace<Scal> ff_flux; // volume flux
    MapCondFaceAdvection<Scal> mf_cond; // face conditions
  } * ctx(sem);

  auto& as = ctx->as;
  auto& mf_cond = ctx->mf_cond;
  auto& fc_src = ctx->fc_src;
  auto& ff_flux = ctx->ff_flux;

  const Scal tmax = 1.;
  const Vect vel(0.5, 0.3, 0.1);
  const double cfl = 0.5;

  if (sem("init")) {
    fc_src.Reinit(m, 0);
    ff_flux.Reinit(m, 0);
    for (auto f : m.Faces()) {
      ff_flux[f] = vel.dot(m.GetSurface(f));
    }
    FieldCell<Scal> fccl(m, 0); // color
    FieldCell<Scal> fcu(m, 0); // initial volume fraction
    for (auto c : m.Cells()) {
      fcu[c] = (m.GetCenter(c).dist(Vect(0.5, 0.5, 0.5)) < 0.2);
    }
    const double dt = cfl * m.GetCellSize()[0] / vel.norm();
    auto p = std::make_shared<typename Vof<M>::Par>();
    {
      Vars vr;
      vr.Double.Set("t0", 0);
      vr.Double.Set("dt", dt);
      vr.Int.Set("max", 10000);
      p->dmp = std::unique_ptr<Dumper>(new Dumper(vr, ""));
    }
    as.reset(new Vof<M>(m, fcu, fccl, mf_cond, &ff_flux, &fc_src, 0., dt, p));
  }
  sem.LoopBegin();
  if (sem.Nested("start")) {
    as->StartStep();
  }
  if (sem.Nested("iter")) {
    as->MakeIteration();
  }
  if (sem.Nested("finish")) {
    as->FinishStep();
  }
  if (sem("checkloop")) {
    if (as->GetTime() >= tmax) {
      sem.LoopBreak();
    }
  }
  if (sem("dump")) {
    m.Dump(&as->GetField(), "u");
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

set int bsx 16
set int bsy 16
set int bsz 16

set int px 1
set int py 1
set int pz 1
)EOF";

  return RunMpiBasic<M, State>(argc, argv, Run, conf);
}
