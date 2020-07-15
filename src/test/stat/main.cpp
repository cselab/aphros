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

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;

void Main(M& m, Vars& var) {
  auto sem = m.GetSem();
  struct {
    FieldNode<Scal> fnl;
    std::unique_ptr<Embed<M>> eb_;
    std::unique_ptr<Stat<M>> stat;
    FieldCell<Scal> vf;
  } * ctx(sem);
  auto& eb_ = ctx->eb_;
  auto& vf = ctx->vf;
  if (sem()) {
    eb_.reset(new Embed<M>(m, 0));
    ctx->fnl = UEmbed<M>::InitEmbed(m, var, false);
  }
  if (sem.Nested()) {
    eb_->Init(ctx->fnl);
  }
  if (sem()) {
    ctx->stat.reset(new Stat<M>(m, eb_.get()));
    auto& stat = *ctx->stat;
    vf.Reinit(m, 0.5);
    stat.Add(
        "vol1", "volume of phase 1", Stat<M>::Reduction::sum,
        [&](IdxCell c, const M& m) { return (1 - vf[c]) * m.GetVolume(c); });
    stat.Add(
        "vol2", "volume of phase 2", Stat<M>::Reduction::sum,
        [&](IdxCell c, const M& m) { return vf[c] * m.GetVolume(c); });
    stat.Add("vol2b", "volume of phase 1", Stat<M>::Reduction::sum, [&]() {
      Scal sum = 0;
      for (auto c : m.Cells()) {
        sum += (1 - vf[c]) * m.GetVolume(c);
      }
      return sum;
    });
    stat.Add(
        "vol2_eb", "volume of phase 1", Stat<M>::Reduction::sum,
        [&](IdxCell c, const Embed<M>& eb) { return vf[c] * eb.GetVolume(c); });
    stat.Add(
        "vf2_max", "max volume fraction 2", Stat<M>::Reduction::max,
        [&](IdxCell c, const M&) { return vf[c]; });
    stat.Add(
        "vf2_max_eb", "max volume fraction 2", Stat<M>::Reduction::max,
        [&](IdxCell c, const Embed<M>&) { return vf[c]; });
    stat.Add(
        "vf2_min", "min volume fraction 2", Stat<M>::Reduction::min,
        [&](IdxCell c, const M&) { return vf[c]; });
  }
  if (sem() && m.IsRoot()) {
    auto& stat = *ctx->stat;
    stat.WriteSummary(std::cout);
    std::cout << std::endl;
    stat.WriteHeader(std::cout);
  }
  for (size_t i = 0; i < 2; ++i) {
    if (sem.Nested()) {
      auto& stat = *ctx->stat;
      stat.Update();
    }
    if (sem() && m.IsRoot()) {
      auto& stat = *ctx->stat;
      stat.WriteValues(std::cout);
    }
    if (sem()) {
      for (auto c : m.CellsM()) {
        vf[c] += c.center[0] * 0.1;
      }
    }
  }
}

int main(int argc, const char** argv) {
  return RunMpiBasicFile<M>(argc, argv, Main);
}
