#undef NDEBUG
#include <cassert>
#include <iostream>
#include <string>
#include <memory>
#include <set>

#include <distr/distrbasic.h>
#include <solver/vofm.h>
#include <debug/isnan.h>
#include <func/init.h>

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;

template <class M>
void CalcMeanCurvatureFlowFlux(
    FieldFace<Scal>& ffv, const GRange<size_t>& layers,
    const Multi<const FieldCell<Scal>*> fcu,
    const Multi<const FieldCell<Scal>*> fccl,
    const Multi<const FieldCell<Vect>*> fcn,
    const Multi<const FieldCell<Scal>*> fck, const M& m) {
  const Scal kClNone = -1;
  Scal vol = 0;
  for (auto c : m.Cells()) {
    for (auto l : layers) {
      vol += (*fcu[l])[c] * m.GetVolume(c);
    }
  }
  vol /= m.GetGlobalLength().prod();
  const Scal vol0 = 0.3;
  ffv.Reinit(m, 0);
  for (auto f : m.Faces()) {
    IdxCell cm = m.GetCell(f, 0);
    IdxCell cp = m.GetCell(f, 1);
    std::set<Scal> s;
    for (auto i : layers) {
      Scal clm = (*fccl[i])[cm];
      Scal clp = (*fccl[i])[cp];
      if (clm != kClNone) s.insert(clm);
      if (clp != kClNone) s.insert(clp);
    }
    for (auto cl : s) {
      Scal um = 0;
      Scal up = 0;
      Vect nm(0);
      Vect np(0);
      Scal km = 0;
      Scal kp = 0;
      for (auto i : layers) {
        if ((*fccl[i])[cm] == cl) {
          um = (*fcu[i])[cm];
          nm = (*fcn[i])[cm];
          km = (*fck[i])[cm];
        }
        if ((*fccl[i])[cp] == cl) {
          up = (*fcu[i])[cp];
          np = (*fcn[i])[cp];
          kp = (*fck[i])[cp];
        }
      }
      const Scal k = (std::abs(um - 0.5) < std::abs(up - 0.5) ? km : kp);
      const Vect n = (std::abs(um - 0.5) < std::abs(up - 0.5) ? nm : np);
      const Scal hi = m.GetArea(f) / m.GetVolume(cp);
      const Scal v = m.GetSurface(f).dot(n * ((vol0 - vol) / vol0)) * hi * k;
      if (!IsNan(v)) {
        ffv[f] = v;
      }
    }
  }
}

void Run(M& m, Vars& var) {
  auto sem = m.GetSem();

  struct {
    std::unique_ptr<Vofm<M>> as; // advection solver
    FieldCell<Scal> fcu; // initial volume fraction
    FieldCell<Scal> fcs; // volume source
    FieldFace<Scal> ffv; // volume flux
    Multi<FieldCell<Scal>> fck; // curvature
    MapCondFaceAdvection<Scal> mf_cond; // face conditions
    GRange<size_t> layers;
  } * ctx(sem);
  auto& fcs = ctx->fcs;
  auto& ffv = ctx->ffv;
  auto& fck = ctx->fck;
  auto& fcu = ctx->fcu;
  auto& layers = ctx->layers;

  auto& as = ctx->as;
  const Scal tmax = 0.1;
  const Scal cfl = 0.1;

  if (sem.Nested()) {
    InitVf(fcu, var, m);
  }
  if (sem("init")) {
    fcs.Reinit(m, 0);
    ffv.Reinit(m, 0);
    FieldCell<Scal> fccl(m, 0); // initial color
    const Scal dt = cfl * m.GetCellSize()[0];
    typename Vofm<M>::Par p;
    as.reset(new Vofm<M>(m, fcu, fccl, ctx->mf_cond, &ffv, &fcs, 0., dt, p));
    layers = GRange<size_t>(as->GetNumLayers());
    fck.resize(layers);
    fck.InitAll(FieldCell<Scal>(m, 1));
  }
  sem.LoopBegin();
  if (sem("flux")) {
    CalcMeanCurvatureFlowFlux(
        ffv, layers, as->GetFieldM(), as->GetColor(), as->GetNormal(), fck, m);
  }
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
set int bsz 1

set int dim 2

set int px 1
set int py 1
set int pz 1

set string init_vf circlels
set vect circle_c 0.5 0.5 0.5
set double circle_r 0.2
)EOF";

  return RunMpiBasic<M>(argc, argv, Run, conf);
}
