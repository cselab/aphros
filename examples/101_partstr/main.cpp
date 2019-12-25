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
    FieldCell<Scal> fccl; // color
    FieldCell<bool> fci; // interface mask
    FieldCell<Scal> fck; // curvature
    typename PSM::Par par;
  } * ctx(sem);
  auto& fca = ctx->fca;
  auto& fcn = ctx->fcn;
  auto& fccl = ctx->fccl;
  auto& fci = ctx->fci;
  auto& fck = ctx->fck;
  auto& par = ctx->par;

  if (sem.Nested()) {
    GRange<size_t> layers(1);
    UCurv<M>::CalcCurvPart(layers, &fca, &fcn, &fci, &fccl, par, &fck, m);
  }
  if (sem("dump")) {
    m.Dump(&fck, "k");
  }
  if (sem()) {
  }
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
