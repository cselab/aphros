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
using Type = typename EB::Type;

void Run(M& m, Vars&) {
  auto sem = m.GetSem("Run");
  struct {
    std::unique_ptr<EB> eb;
    FieldCell<Scal> fcu;
    FieldCell<Vect> fcg;
    FieldEmbed<Scal> feu;
    size_t frame = 0;
  } * ctx(sem);
  auto& fcu = ctx->fcu;
  auto& fcg = ctx->fcg;
  auto& feu = ctx->feu;
  auto& frame = ctx->frame;

  if (sem("init")) {
    FieldNode<Scal> fnl(m, 0);
    for (auto n : m.AllNodes()) {
      const Vect x = m.GetNode(n);
      const Vect xc(0.5, 0.5, 0.5);
      const Vect s(1., 1., 1.);
      auto rot = [](Vect xx) {
        const Scal a = 3.14159265358979323 * 0;
        const Scal as = std::sin(a);
        const Scal ac = std::cos(a);
        const Scal x = xx[0];
        const Scal y = xx[1];
        const Scal z = xx[2];
        return Vect(x * ac - y * as, x * as + y * ac, z);
      };
      fnl[n] = (rot(x - xc) / s).norm() - 0.35;
      fnl[n] *= -1;
    }
    ctx->eb.reset(new EB(m, fnl));
    fcu.Reinit(m, 0);
  }
  if (sem.Nested("dumppoly")) {
    auto& eb = *ctx->eb;
    eb.DumpPoly();
  }
  const size_t maxt = 100;
  const size_t nfr = 100;
  for (size_t t = 0; t < maxt; ++t) {
    auto& eb = *ctx->eb;
    if (sem("step")) {
      feu = eb.Interpolate(fcu, 0, 1);

      /*
      const Vect g(0.3, 0.2, 0.1);
      for (auto c : m.AllCells()) {
        feu[c] = eb.GetFaceCenter(c).dot(g);
      }
      for (auto f : m.AllFaces()) {
        feu[f] = eb.GetFaceCenter(f).dot(g);
      }
      */

      for (auto c : m.AllCells()) {
        fcu[c] = Vect(0.5).dist(eb.GetCellCenter(c));
      }
      feu = eb.Gradient(fcu, 0, 0.35);

      fcu = eb.Interpolate(feu);
      fcg = eb.Gradient(feu);
      m.Comm(&fcu);
    }
    if (t % std::max<size_t>(1, maxt / nfr) != 0) {
      continue;
    }
    if (sem("dumpcsv")) {
      auto& eb = *ctx->eb;
      std::ofstream out(GetDumpName("eb", ".csv", frame));
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
      std::ofstream out(GetDumpName("ebf", ".csv", frame));
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
    if (sem("dump")) {
      // FIXME: Dump and Comm in one stage ignores Comm
      m.Dump(&fcu, "u"); 
      m.Dump(&fcg, 0, "ux"); 
      m.Dump(&fcg, 1, "uy"); 
      m.Dump(&fcg, 2, "uz"); 
      ++frame;
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

set int bsx 16
set int bsy 16
set int bsz 16
)EOF";
  return RunMpiBasic<M>(argc, argv, Run, conf);
}
