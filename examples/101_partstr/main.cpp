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
  auto sem = m.GetSem();

  struct {
    FieldFace<Scal> fca; // plane constant
    FieldFace<Vect> fcn; // normal
    FieldFace<Scal> fccl; // color
    FieldFace<bool> fci; // interface mask
    FieldFace<Scal> fck; // curvature
    typename PartstrMeshM<M>::Par par;
    GRange<size_t> layers(1);
  } * ctx(sem);
  auto& fca = ctx->fca;
  auto& fcn = ctx->fcn;
  auto& fccl = ctx->fccl;
  auto& fci = ctx->fci;
  auto& par = ctx->par;

  if (sem.Nested()) {
    UCurv<M>::CalcCurvPart(layers, &fca, &fcn, &fci, &fccl, par, fck, m);
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
