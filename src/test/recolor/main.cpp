#undef NDEBUG
#include <iostream>

#include "distr/distrbasic.h"
#include "func/init_u.h"
#include "solver/solver.h"
#include "util/vof.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;

struct State {
  FieldCell<Scal> fcu;
  FieldCell<Scal> fccl;
  solver::UVof<M> uvof;
  MapFace<std::shared_ptr<solver::CondFace>> mfc;
};

void Run(M& m, State& s, Vars& var) {
  (void) var;
  auto sem = m.GetSem();
  auto& fcu = s.fcu;
  auto& fccl = s.fccl;
  auto& uvof = s.uvof;
  auto& mfc = s.mfc;
  constexpr Scal kClNone = -1;
  GRange<size_t> layers{0, 1};
  if (sem()) {
    fcu.Reinit(m);
    auto init = CreateInitU<M>(var, m.IsRoot());
    init(fcu, m);
    fccl.Reinit(m, kClNone);
  }
  for (size_t i = 0; i < 10; ++i) {
    if (sem()) {
      for (auto c : m.Cells()) {
        fccl[c] = (fcu[c] > 0 ? 1. : kClNone);
      }
      m.Comm(&fcu);
      m.Comm(&fccl);
    }
    if (sem.Nested()) {
      uvof.Recolor(layers, &fcu, &fccl, -1, Vect(0), 1e10, mfc, true, true, m);
    }
    if (sem()) {
      for (auto c : m.Cells()) {
        auto& cl = fccl[c];
        if (cl != kClNone) {
          if (cl == 0) {
            // nop
          } else {
            cl = 1. + std::sin(1234567 * cl);
          }
        }
      }
      m.Dump(&fcu, "u");
      m.Dump(&fccl, "cl");
    }
    if (sem.Nested()) {
      solver::Smoothen(fcu, mfc, m, 1);
    }
  }
  if (sem()) {}
}

int main(int argc, const char** argv) {
  std::string conf = R"EOF(
set int bx 4
set int by 4
set int bz 2

set int bsx 8
set int bsy 8
set int bsz 8

set string init_vf list
set int list_ls 1
set string list_path b.dat
  )EOF";

  return RunMpiBasic<M, State>(argc, argv, Run, conf);
}

