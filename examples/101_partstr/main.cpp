#undef NDEBUG
#include <cassert>
#include <iostream>
#include <memory>
#include <string>

#include "distr/distrbasic.h"
#include "solver/curv.h"

using namespace solver;

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using namespace solver;

struct State {};

void Run(M& m, State&, Vars&) {
  using PSM = PartStrMeshM<M>;
  auto sem = m.GetSem();

  struct {
    FieldCell<Scal> fca; // plane constant
    FieldCell<Vect> fcn; // normal
    FieldCell<bool> fci; // interface mask
    FieldCell<Scal> fck; // curvature
    typename PSM::Par par;
    std::unique_ptr<PSM> psm;
  } * ctx(sem);
  auto& fca = ctx->fca;
  auto& fcn = ctx->fcn;
  auto& fci = ctx->fci;
  auto& fck = ctx->fck;
  auto& par = ctx->par;
  auto& psm = ctx->psm;

  if (sem("init")) {
    auto ps = std::make_shared<typename PartStr<Scal>::Par>();
    par.ps = ps;
    fca.Reinit(m, 0);
    fcn.Reinit(m, Vect(0.5));
    fci.Reinit(m, true);
  }
  if (sem.Nested()) {
    psm = UCurv<M>::CalcCurvPart(
        GRange<size_t>(1), &fca, &fcn, &fci, nullptr, par, &fck, m);
  }
  if (sem.Nested()) {
    psm->DumpParticles(&fca, &fcn, 0, 0);
  }
  if (sem("dump")) {
    m.Dump(&fck, "k");
  }
  if (sem()) {}
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
