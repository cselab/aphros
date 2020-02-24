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
    MapEmbed<BCondFluid<Vect>> mebc;
    MapEmbed<size_t> megroup;
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
    const std::string filename = "bc.dat";
    std::ifstream fin(filename);
    if (!fin.good()) {
      throw std::runtime_error(
          "Can't open list of boundary conditions '" + filename + "'");
    }
    std::vector<std::string> vdesc;
    std::tie(ctx->mebc, ctx->megroup, vdesc) = UInitEmbedBc<M>::Parse(fin, eb);
    fin.close();
    std::ofstream fdesc("bc_groups.dat");
    for (size_t i = 0; i < vdesc.size(); ++i) {
      fdesc << i << " " << vdesc[i] << std::endl;
    }
    std::ofstream out("bc.csv");
    out << "x,y,z,type,group\n";
    auto write = [&out](const BCondFluid<Vect>& bc, Vect x, size_t g) {
      out << x[0] << "," << x[1] << "," << x[2];
      out << "," << int(bc.type) << "," << g << "\n";
    };
    for (auto& p : ctx->mebc.GetMapCell()) {
      write(p.second, eb.GetFaceCenter(p.first), ctx->megroup[p.first]);
    }
    for (auto& p : ctx->mebc.GetMapFace()) {
      write(p.second, eb.GetFaceCenter(p.first), ctx->megroup[p.first]);
    }
  }
  if (sem.Nested("bcdump")) {
    UInitEmbedBc<M>::DumpPoly("bc.vtk", ctx->megroup, *ctx->eb, m);
  }
  const size_t maxt = 10;
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

      const auto feunc = UEB::GradientBilinear(fcu, mfc, bc, bcu, eb);
      if (compact) {
        feun = feunc;
      } else {
        const FieldCell<Vect> fcg = UEB::GradientGauss(feu, eb);
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
set int bx 2
set int by 2
set int bz 2

set int bsx 16
set int bsy 16
set int bsz 16

set string eb_init box
set vect eb_box_c 0.5 0.5 0.5
#set vect eb_box_r 0.249 0.249 0.249
set vect eb_box_r 10 0.249 10

set string eb_init sphere
set vect eb_sphere_c 0.3 0.5 0.5
set vect eb_sphere_r 0.5 0.349 0.249
set double eb_sphere_angle 0
set int eb_init_inverse 1
)EOF";
  return RunMpiBasic<M>(argc, argv, Run, conf);
}
