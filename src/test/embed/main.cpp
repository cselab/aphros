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
    FieldEmbed<Scal> feu;
  } * ctx(sem);
  auto& fcu = ctx->fcu;
  auto& fcum = ctx->fcum;
  auto& feu = ctx->feu;

  if (sem("init")) {
    FieldNode<Scal> fnl(m, 0);
    for (auto n : m.AllNodes()) {
      const Vect x = m.GetNode(n);
      const Vect xc(0.5, 0.5, 0.5);
      fnl[n] = (x - xc).norminf() - 0.2;
    }
    ctx->eb.reset(new EB(m, fnl));
    fcu.Reinit(m, 0);
  }
  if (sem.Nested("dumppoly")) {
    auto& eb = *ctx->eb;
    eb.DumpPoly();
  }
  for (size_t t = 0; t < 50; ++t) {
    auto& eb = *ctx->eb;
    if (sem("step")) {
      feu = eb.Interpolate(fcu, 0, 1);
      fcu = eb.Interpolate(feu);
      m.Comm(&fcu);
      m.Dump(&fcu, "u");
    }
    if (sem("dumpcsv")) {
      auto& eb = *ctx->eb;
      using Type = typename EB::Type;
      std::ofstream out(GetDumpName("eb", ".csv", t));
      out << "x,y,z,type,u\n";
      for (auto c : m.Cells()) {
        if (eb.GetType(c) == Type::regular || eb.GetType(c) == Type::cut) {
          auto x = eb.GetCellCenter(c);
          out << x[0] << "," << x[1] << "," << x[2];
          out << "," << size_t(eb.GetType(c));
          out << "," << fcu[c] << "\n";
        }
      }
    }
    if (sem("dumpcsvface")) {
      auto& eb = *ctx->eb;
      using Type = typename EB::Type;
      std::ofstream out(GetDumpName("ebf", ".csv", t));
      out << "x,y,z,face,type,u\n";
      for (auto c : m.Cells()) {
        if (eb.GetType(c) == Type::cut) {
          auto x = eb.GetFaceCenter(c);
          out << x[0] << "," << x[1] << "," << x[2];
          out << "," << 0;
          out << "," << size_t(eb.GetType(c));
          out << "," << feu[c] << "\n";
        }
      }
      for (auto f : m.Faces()) {
        if (eb.GetType(f) == Type::regular || eb.GetType(f) == Type::cut) {
          auto x = eb.GetFaceCenter(f);
          out << x[0] << "," << x[1] << "," << x[2];
          out << "," << 1;
          out << "," << size_t(eb.GetType(f));
          out << "," << feu[f] << "\n";
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
