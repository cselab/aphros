// Created by Petr Karnakov on 15.07.2020
// Copyright 2020 ETH Zurich

#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>

#include "distr/distrbasic.h"
#include "geom/mesh.h"
#include "parse/vars.h"
#include "solver/approx_eb.h"
#include "solver/embed.h"
#include "util/stat.h"

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;

void Main(M& m, Vars& var) {
  auto sem = m.GetSem();
  struct {
    FieldNode<Scal> fnl;
    std::unique_ptr<Embed<M>> eb_;
    std::unique_ptr<Stat<M>> stat;
    FieldCell<Scal> vf;
    std::ofstream fout;
  } * ctx(sem);
  auto& eb_ = ctx->eb_;
  auto& vf = ctx->vf;
  if (sem()) {
    eb_.reset(new Embed<M>(m, 0));
  }
  if (sem.Nested("levelset")) {
    UEmbed<M>::InitLevelSet(ctx->fnl, m, var, false);
  }
  if (sem.Nested()) {
    eb_->Init(ctx->fnl);
  }
  if (sem()) {
    ctx->stat.reset(new Stat<M>(m, eb_.get()));
    auto& stat = *ctx->stat;
    vf.Reinit(m, 0.5);
    stat.AddSum("vol", "volume", [&](IdxCell c, const M&) { //
      return m.GetVolume(c);
    });
    stat.AddSumHidden("xvfvol2", "x*vf*volume", [&](IdxCell c, const M&) { //
      return m.GetCenter(c) * m.GetVolume(c) * vf[c];
    });
    stat.AddSum("vol1", "volume of phase 1", [&](IdxCell c, const M&) {
      return (1 - vf[c]) * m.GetVolume(c);
    });
    stat.AddSum("vol2", "volume of phase 2", [&](IdxCell c, const M&) {
      return vf[c] * m.GetVolume(c);
    });
    stat.AddSum("vol2b", "volume of phase 1", [&]() {
      Scal sum = 0;
      for (auto c : m.Cells()) {
        sum += (1 - vf[c]) * m.GetVolume(c);
      }
      return sum;
    });
    stat.AddSum(
        "vol2_eb", "volume of phase 1",
        [&](IdxCell c, const Embed<M>& eb) { return vf[c] * eb.GetVolume(c); });
    stat.AddMax("vf2_max", "max volume fraction 2", [&](IdxCell c, const M&) {
      return vf[c];
    });
    stat.AddMax(
        "vf2_max_eb", "max volume fraction 2",
        [&](IdxCell c, const Embed<M>&) { return vf[c]; });
    stat.AddMin("vf2_min", "min volume fraction 2", [&](IdxCell c, const M&) {
      return vf[c];
    });
    stat.AddDerived("vol_copy", "copy of volume", [](const Stat<M>& s) {
      return s["vol"];
    });
    stat.AddDerived("c2", "centeroid of phase 2", [](const Stat<M>& s) {
      return s.vect["xvfvol2"] / s["vol2"];
    });
  }
  if (sem() && m.IsRoot()) {
    auto& stat = *ctx->stat;
    stat.WriteSummary(std::cout);
    std::cout << std::endl;
    ctx->fout.open("stat.dat");
    stat.WriteHeader(ctx->fout);
  }
  for (size_t i = 0; i < 2; ++i) {
    if (sem.Nested()) {
      auto& stat = *ctx->stat;
      stat.Update();
    }
    if (sem() && m.IsRoot()) {
      auto& stat = *ctx->stat;
      stat.WriteValues(ctx->fout);
    }
    if (sem()) {
      for (auto c : m.CellsM()) {
        vf[c] += Vect(0.1, 0.2, 0.3).dot(c.center);
      }
    }
  }
}

int main(int argc, const char** argv) {
  MpiWrapper mpi(&argc, &argv);
  return RunMpiBasicFile<M>(mpi, Main);
}
