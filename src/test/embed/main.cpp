#undef NDEBUG
#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <utility>

#include "distr/distrbasic.h"
#include "solver/reconst.h"
#include "solver/embed.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using R = Reconst<Scal>;
using EB = Embed<M>;

void Run(M& m, Vars&) {
  auto sem = m.GetSem("Run");
  struct {
    std::unique_ptr<EB> eb;
    FieldCell<Scal> fcu;
    FieldCell<Scal> fcum;
    FieldCell<Scal> fct;
    FieldEmbed<Scal> feu;
  } * ctx(sem);
  auto& fct = ctx->fct;
  auto& fcu = ctx->fcu;
  auto& fcum = ctx->fcum;
  auto& feu = ctx->feu;

  if (sem("init")) {
    FieldNode<Scal> fnf(m, 0);
    for (auto n : m.Nodes()) {
      auto x = m.GetNode(n);
      fnf[n] = 1.01 - Vect(x[0], x[1], x[2]).dot(Vect(1., 1., 0.));
    }
    ctx->eb.reset(new EB(m, fnf));
    fcu.Reinit(m, 0);
    fct.Reinit(m, 0);
    for (auto c : m.Cells()) {
      fct[c] = size_t(ctx->eb->GetType(c));
    }
  }
  if (sem.Nested("dumppoly")) {
    auto& eb = *ctx->eb;
    eb.DumpPoly();
  }
  for (size_t t = 0; t < 50; ++t) {
    auto& eb = *ctx->eb;
    if (sem("step")) {
      using Type = typename EB::Type;
      const Scal a = 1.; // value on boundary
      const Vect vel(1., 0., 0.); // advection velocity
      const Scal dt = 0.1 * m.GetCellSize()[0] / vel.norm();

      feu = eb.Interpolate(fcu, 0, 1);
      fcu = eb.Interpolate(feu);
      m.Comm(&fcu);
      m.Dump(&fcu, "u");
    }
    if (sem("dumpcsv")) {
      auto& eb = *ctx->eb;
      using Type = typename EB::Type;
      std::ofstream out(GetDumpName("eb", ".csv", t));
      out << "x,y,z,u\n";
      for (auto c : m.Cells()) {
        if (eb.GetType(c) == Type::cut) {
          auto x = eb.GetCellCenter(c);
          out << x[0] << "," << x[1] << "," << x[2] << "," << fcu[c] << "\n";
        }
      }
    }
    if (sem("dumpcsvface")) {
      auto& eb = *ctx->eb;
      using Type = typename EB::Type;
      std::ofstream out(GetDumpName("ebf", ".csv", t));
      out << "x,y,z,u\n";
      for (auto c : m.Cells()) {
        if (eb.GetType(c) == Type::cut) {
          auto x = eb.GetFaceCenter(c);
          out << x[0] << "," << x[1] << "," << x[2] << "," << feu[c] << "\n";
        }
      }
      for (auto f : m.Faces()) {
        if (eb.GetType(f) == Type::cut) {
          auto x = eb.GetFaceCenter(f);
          out << x[0] << "," << x[1] << "," << x[2] << "," << feu[f] << "\n";
        }
      }
    }
  }
  if (sem()) {
  }
}

int main(int argc, const char** argv) {
  std::string conf = R"EOF(
set int bx 1
set int by 1
set int bz 1

set int bsx 8
set int bsy 8
set int bsz 8
)EOF";
  return RunMpiBasic<M>(argc, argv, Run, conf);
}
