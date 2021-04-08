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
#include "util/hydro.h"

using M = MeshCartesian<double, 3>;
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
    FieldNode<Scal> fnl;
    size_t frame = 0;
    MapEmbed<BCondFluid<Vect>> mebc;
    MapEmbed<size_t> me_group;
    MapEmbed<Scal> me_contang;
  } * ctx(sem);
  auto& fcu = ctx->fcu;
  auto& feu = ctx->feu;
  auto& frame = ctx->frame;

  if (sem("ctor")) {
    ctx->eb.reset(new EB(m));
    fcu.Reinit(m, 0);
    feu.Reinit(m, 0);
  }
  if (sem.Nested("levelset")) {
    UEB::InitLevelSet(ctx->fnl, m, var, m.IsRoot());
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
    fassert(
        fin.good(),
        "Can't open list of boundary conditions '" + filename + "'");
    std::vector<std::string> vdesc;
    MapEmbed<size_t> me_nci;
    std::tie(ctx->me_group, me_nci, vdesc) =
        UInitEmbedBc<M>::ParseGroups(fin, eb, FieldCell<bool>());
    fin.close();

    ctx->me_group.LoopPairs([&](auto p) {
      const auto cf = p.first;
      ctx->me_contang[cf] = 0;
    });

    ctx->mebc =
        UInitEmbedBc<M>::GetBCondFromGroups(ctx->me_group, me_nci, vdesc, eb);

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
      write(p.second, eb.GetFaceCenter(p.first), ctx->me_group[p.first]);
    }
    for (auto& p : ctx->mebc.GetMapFace()) {
      write(p.second, eb.GetFaceCenter(p.first), ctx->me_group[p.first]);
    }
  }
  if (sem.Nested("bcdump")) {
    DumpBcPoly("bc.vtk", ctx->me_group, ctx->me_contang, *ctx->eb, m);
  }
  const size_t maxt = 10;
  for (size_t t = 0; t < maxt; ++t) {
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

  MpiWrapper mpi(&argc, &argv);
  return RunMpiBasicString<M>(mpi, Run, conf);
}
