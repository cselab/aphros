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
  Scal sum;
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
    fccl.Reinit(m, 1);
  }
  for (size_t i = 0; i < 10; ++i) {
    if (sem()) {
      for (auto c : m.Cells()) {
        if (fcu[c] == 0) {
          fccl[c] = kClNone;
        } else {
          if (fccl[c] == kClNone) {
            for (auto q : m.Nci(c)) {
              auto cn = m.GetCell(c, q);
              if (fccl[cn] != kClNone) {
                fccl[c] = fccl[cn];
                break;
              }
            }
          }
        }
      }
      m.Comm(&fcu);
      m.Comm(&fccl);
    }
    if (sem.Nested()) {
      uvof.Recolor(layers, &fcu, &fccl, -1, Vect(0), 1e10, mfc, true, false, m);
    }
    if (sem()) {
      s.sum = 0;
      for (auto c : m.Cells()) {
        s.sum += fcu[c] * m.GetVolume(c);
      }
      m.Reduce(&s.sum, "sum");
    }
    if (sem()) {
      /*
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
      */
      if (m.IsRoot()) {
        std::cout << "sum=" << s.sum << std::endl;
      }
      m.Dump(&fcu, "u");
      m.Dump(&fccl, "cl");
    }
    if (sem.Nested()) {
      solver::Smoothen(fcu, mfc, m, 1);
    }
  }
  if (sem()) {} // empty stage for dump
}

int main(int argc, const char** argv) {
  std::string conf = R"EOF(
set int bx 2
set int by 2
set int bz 1

set int bsx 16
set int bsy 16
set int bsz 16

set string init_vf list
set int list_ls 1
set string list_path b.dat
  )EOF";

  return RunMpiBasic<M, State>(argc, argv, Run, conf);
}

