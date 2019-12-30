#undef NDEBUG
#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <utility>

#include "distr/distrbasic.h"
#include "solver/embed.h"
#include "solver/reconst.h"

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
    FieldEmbed<Scal> feu;
    size_t frame = 0;
  } * ctx(sem);
  auto& fcu = ctx->fcu;
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
  const size_t maxt = 1000;
  const size_t nfr = 100;
  for (size_t t = 0; t < maxt; ++t) {
    auto& eb = *ctx->eb;
    if (sem("step")) {
      const bool compact = true;
      //const Scal dt = compact ? 6e-6 : 4e-5;
      const Scal dt = compact ? 1e-4 : 1e-3; // with redistr
      const Scal bcu = 1;

      feu = eb.Interpolate(fcu, 0, bcu);
      FieldEmbed<Scal> feun(m);

      if (compact) {
        feun = eb.Gradient(fcu, 0, bcu); // normal gradient
      } else {
        const FieldCell<Vect> fcg = eb.Gradient(feu);
        const FieldEmbed<Vect> feg = eb.Interpolate(fcg, 1, Vect(0));
        for (auto f : m.Faces()) {
          if (eb.GetType(f) != Type::excluded) {
            feun[f] = feg[f].dot(eb.GetNormal(f));
          }
        }
        for (auto c : m.Cells()) {
          if (eb.GetType(c) == Type::cut) {
            feun[c] = feg[c].dot(eb.GetNormal(c));
          }
        }
      }

      for (auto c : m.Cells()) {
        if (eb.GetType(c) == Type::excluded) {
          continue;
        }
        Scal s = 0;
        if (eb.GetType(c) == Type::cut) {
          s += feun[c] * eb.GetArea(c);
        }
        for (auto q : m.Nci(c)) {
          auto f = m.GetFace(c, q);
          s += feun[f] * eb.GetArea(f) * m.GetOutwardFactor(c, q);
        }
        const Scal du = s * dt / eb.GetVolume(c);
        if (eb.GetType(c) == Type::cut) {
          fcu[c] += du * eb.GetRedistr(c);
          for (auto p : eb.GetRedistrList(c)) {
            fcu[p.first] += du * p.second;
          }
        } else {
          fcu[c] += du;
        }
      }

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
