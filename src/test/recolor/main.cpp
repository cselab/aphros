// Created by Petr Karnakov on 14.10.2019
// Copyright 2019 ETH Zurich

#undef NDEBUG
#include <iostream>
#include <memory>

#include "distr/distrbasic.h"
#include "dump/dump.h"
#include "func/init_u.h"
#include "solver/approx_eb.h"
#include "solver/solver.h"
#include "solver/trackerm.h"
#include "util/vof.h"

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using TRM = Trackerm<M>;

void Dump(M& m, const FieldCell<Scal>& fc, std::string name, int i) {
  size_t id = m.GetId();
  id = id % 1024 + (id / 1024) * 2;
  std::string sid = std::to_string(id);
  std::string op = name + "_" + sid + "_" + std::to_string(i) + ".dat";
  dump::DumpFieldPlain(fc, m.GetIndexCells(), m.GetAllBlockCells(), op);
}

void Run(M& m, Vars& var) {
  auto sem = m.GetSem();
  struct {
    FieldCell<Scal> fcu;
    FieldCell<Scal> fccl;
    FieldCell<Scal> fcclm;
    UVof<M> uvof;
    MapEmbed<BCond<Scal>> mfc;
    std::unique_ptr<TRM> trm;
    FieldCell<Scal> fcim;
  } * ctx(sem);

  auto& fcu = ctx->fcu;
  auto& fccl = ctx->fccl;
  auto& fcclm = ctx->fcclm;
  auto& uvof = ctx->uvof;
  auto& mfc = ctx->mfc;
  constexpr Scal kClNone = -1;
  GRange<size_t> layers(1);

  auto Recolor = [&]() {
    if (sem.Nested()) {
      ctx->trm->Update(&fccl, &fcclm);
    }
    if (sem.Nested()) {
      const bool unionfind = true;
      const bool reduce = true;
      const bool grid = true;
      uvof.Recolor(
          layers, &fcu, &fccl, &fccl, -1, Vect(0), 1e10, mfc, true, unionfind,
          reduce, grid, m);
      fcclm = fccl;
    }
  };

  if (sem.Nested()) {
    InitVf(fcu, var, m, true);
  }
  if (sem()) {
    fccl.Reinit(m);
    ctx->fcim.Reinit(m);
    for (auto c : m.Cells()) {
      fccl[c] = (fcu[c] == 0 ? kClNone : 1);
    }
    fcclm = fccl;
    m.Comm(&fcu);
    m.Comm(&fccl);
    ctx->trm = std::unique_ptr<TRM>(new TRM(m, layers));
  }
  Recolor();
  for (size_t i = 0; i < 5; ++i) {
    if (sem()) {
      m.Dump(&fcu, "u");
      m.Dump(&fccl, "cl");
      for (auto c : m.Cells()) {
        ctx->fcim[c] = ctx->trm->GetImage(0, c)[1];
      }
      m.Dump(&ctx->fcim, "im");
    }
    if (sem.Nested()) {
      Smoothen(fcu, mfc, m, 1);
    }
    if (sem()) {
      fcclm = fccl;
      for (auto c : m.Cells()) {
        if (fcu[c] == 0) {
          fccl[c] = kClNone;
        } else {
          if (fcclm[c] == kClNone) {
            for (auto q : m.Nci(c)) {
              auto cn = m.GetCell(c, q);
              if (fcclm[cn] != kClNone) {
                fccl[c] = fcclm[cn];
                break;
              }
            }
          }
        }
      }
      m.Comm(&fccl);
    }
    Recolor();
  }
  if (sem()) {
  } // empty stage for dump
}

int main(int argc, const char** argv) {
  std::string conf = R"EOF(
set int bx 4
set int by 4
set int bz 4

set int bsx 8
set int bsy 8
set int bsz 8

set int openmp 0
set int histogram 0
set int mpi_compress_msg 0

set string init_vf list
set int list_ls 1
set string list_path b.dat
  )EOF";

  MpiWrapper mpi(&argc, &argv);
  return RunMpiBasicString<M>(mpi, Run, conf);
}
