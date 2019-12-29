#undef NDEBUG
#include <cassert>
#include <iostream>
#include <memory>
#include <set>
#include <string>

#include <debug/isnan.h>
#include <distr/distrbasic.h>
#include <func/init.h>
#include <solver/curv.h>
#include <solver/vofm.h>

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
  const Scal kvol = 100.;
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
      const Scal v = m.GetSurface(f).dot(n * (-k + kvol * (vol0 - vol) / vol0));
      if (!IsNan(v)) {
        ffv[f] = v;
      }
      /*
      const Scal k = (std::abs(um - 0.5) < std::abs(up - 0.5) ? km : kp);
      const Scal hi = m.GetArea(f) / m.GetVolume(cp);
      const Scal ga = -(up - um) * hi * (-k + kvol * (vol0 - vol) / vol0);
      if (std::abs(ga) > 0.) {
        ffv[f] += ga * m.GetArea(f);
      }
      */
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
    typename PartStrMeshM<M>::Par psm_par;
    std::unique_ptr<PartStrMeshM<M>> psm;
    Scal dt = 0.01;
    size_t step = 0;
  } * ctx(sem);
  auto& fcs = ctx->fcs;
  auto& ffv = ctx->ffv;
  auto& fck = ctx->fck;
  auto& fcu = ctx->fcu;
  auto& layers = ctx->layers;
  auto& psm_par = ctx->psm_par;
  auto& psm = ctx->psm;
  auto& dt = ctx->dt;
  auto& step = ctx->step;

  auto& as = ctx->as;
  const Scal tmax = 0.1;
  const Scal cfl = 0.5;

  if (sem.Nested()) {
    InitVf(fcu, var, m);
  }
  if (sem("init")) {
    fcs.Reinit(m, 0);
    ffv.Reinit(m, 0);
    FieldCell<Scal> fccl(m, 0); // initial color
    typename Vofm<M>::Par p;
    as.reset(new Vofm<M>(m, fcu, fccl, ctx->mf_cond, &ffv, &fcs, 0., dt, p));
    layers = GRange<size_t>(as->GetNumLayers());
    fck.resize(layers);
    fck.InitAll(FieldCell<Scal>(m, 1));
    psm_par.dump_fr = 1;
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
    ++step;
  }
  if (sem.Nested()) {
    psm = UCurv<M>::CalcCurvPart(
        layers, as->GetAlpha(), as->GetNormal(), as->GetMask(), as->GetColor(),
        psm_par, fck, m);
  }
  if (sem("flux")) {
    CalcMeanCurvatureFlowFlux(
        ffv, layers, as->GetFieldM(), as->GetColor(), as->GetNormal(), fck, m);
    Scal maxv = 0;
    for (auto f : m.Faces()) {
      maxv = std::max(maxv, std::abs(ffv[f]));
    }
    dt = cfl * m.GetCellSize().prod() / maxv;
    as->SetTimeStep(dt);
    if (m.IsRoot()) {
      std::cout << "dt=" << dt << std::endl;
    }
  }
  if (sem.Nested()) {
    psm->DumpParticles(as->GetAlpha(), as->GetNormal(), step, as->GetTime());
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
