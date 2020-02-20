// Created by Petr Karnakov on 13.05.2019
// Copyright 2019 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <utility>

#include "distr/distrbasic.h"
#include "func/init_bc.h"
#include "solver/approx_eb.h"
#include "solver/embed.h"
#include "solver/reconst.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using R = Reconst<Scal>;
using EB = Embed<M>;
using UEB = UEmbed<M>;
using Type = typename EB::Type;

void Run(M& m, Vars& var) {
  auto sem = m.GetSem("Run");
  struct {
    std::unique_ptr<EB> eb;
    FieldCell<Scal> fcu;
    FieldEmbed<Scal> feu;
    MapCondFace mfc;
    FieldNode<Scal> fnl;
    size_t frame = 0;
  } * ctx(sem);
  auto& fcu = ctx->fcu;
  auto& feu = ctx->feu;
  auto& frame = ctx->frame;
  auto& mfc = ctx->mfc;

  if (sem("ctor")) {
    ctx->eb.reset(new EB(m));
    ctx->fnl = UEB::InitEmbed(m, var, m.IsRoot());
    fcu.Reinit(m, 0);
  }
  if (sem.Nested("init")) {
    ctx->eb->Init(ctx->fnl);
  }
  if (sem.Nested("dumppoly")) {
    auto& eb = *ctx->eb;
    eb.DumpPoly();
  }
  if (sem("bc")) {
    using namespace fluid_condition;
    auto& eb = *ctx->eb;
    MapEmbed<UniquePtr<CondFaceFluid>> mec;
    const std::string fn = "bc.dat";
    UInitEmbedBc<M>::Parse(fn, eb, mec);
    std::ofstream out("bc.csv");
    out << "x,y,z,bc\n";
    for (auto c : eb.CFaces()) {
      auto& b = mec[c];
      Scal condf = -1;
      if (!b.Get()) {
        condf = 0;
      } else if (b.Get<NoSlipWall<M>>()) {
        condf = 1;
      } else if (b.Get<SlipWall<M>>()) {
        condf = 2;
      } else if (b.Get<Inlet<M>>()) {
        condf = 3;
      } else if (b.Get<InletFlux<M>>()) {
        condf = 4;
      } else if (b.Get<Outlet<M>>()) {
        condf = 5;
      } else if (b.Get<Symm<M>>()) {
        condf = 6;
      }
      auto x = eb.GetCellCenter(c);
      out << x[0] << "," << x[1] << "," << x[2] << "," << condf << "\n";
    }
  }
  const size_t maxt = 100;
  const size_t nfr = 10;
  for (size_t t = 0; t < maxt; ++t) {
    auto& eb = *ctx->eb;
    if (sem("step")) {
      const bool compact = false;
      // const Scal dt = compact ? 6e-6 : 4e-5;
      const Scal dt = compact ? 1e-4 : 5e-4; // with redistr
      // const Scal dt = compact ? 5e-4 : 5e-4; // with redistr, bc=1
      const size_t bc = 0;
      const Scal bcu = 1;

      feu = UEB::Interpolate(fcu, mfc, bc, bcu, eb);
      FieldEmbed<Scal> feun(m);

      const auto feunc = UEB::Gradient(fcu, mfc, bc, bcu, eb);
      if (compact) {
        feun = feunc;
      } else {
        const FieldCell<Vect> fcg = UEB::Gradient(feu, eb);
        const FieldEmbed<Vect> feg =
            UEB::Interpolate(fcg, MapCondFace(), 1, Vect(0), eb);
        for (auto f : eb.Faces()) {
          feun[f] = feg[f].dot(eb.GetNormal(f));
        }
        for (auto c : eb.CFaces()) {
          feun[c] = feg[c].dot(eb.GetNormal(c));
          const Scal a = eb.GetRedistr(c);
          feun[c] = feunc[c] * a + feun[c] * (1 - a);
        }
      }

      for (auto c : eb.Cells()) {
        Scal s = 0;
        if (eb.GetType(c) == Type::cut) {
          s += feun[c] * eb.GetArea(c);
        }
        for (auto q : eb.Nci(c)) {
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
      for (auto c : eb.Cells()) {
        auto x = eb.GetCellCenter(c);
        out << x[0] << "," << x[1] << "," << x[2];
        out << "," << size_t(eb.GetType(c));
        out << "," << fcu[c] << "\n";
      }
    }
    if (sem("dumpcsvface")) {
      auto& eb = *ctx->eb;
      std::ofstream out(GetDumpName("ebf", ".csv", frame));
      out << "x,y,z,face,type,u\n";
      for (auto c : eb.CFaces()) {
        auto x = eb.GetFaceCenter(c);
        out << x[0] << "," << x[1] << "," << x[2];
        out << "," << 0;
        out << "," << size_t(eb.GetType(c));
        out << "," << feu[c] << "\n";
      }
      for (auto f : eb.Faces()) {
        auto x = eb.GetFaceCenter(f);
        out << x[0] << "," << x[1] << "," << x[2];
        out << "," << 1;
        out << "," << size_t(eb.GetType(f));
        out << "," << feu[f] << "\n";
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

set string eb_init box
set vect eb_box_c 0.5 0.5 0.5
#set vect eb_box_r 0.249 0.249 0.249
set vect eb_box_r 10 0.249 10

set string eb_init sphere
set vect eb_sphere_c 0.5 0.5 0.5
set vect eb_sphere_r 0.249 0.249 0.249
set double eb_sphere_angle 0
set int eb_init_inverse 0
)EOF";
  return RunMpiBasic<M>(argc, argv, Run, conf);
}
