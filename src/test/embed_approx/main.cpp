// Created by Petr Karnakov on 05.07.2020
// Copyright 2020 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

#include <distr/distrbasic.h>
#include "dump/dump.h"
#include "kernel/kernelmeshpar.h"
#include "parse/vars.h"
#include "solver/approx.h"
#include "solver/approx_eb.h"
#include "solver/approx_eb.ipp"
#include "solver/embed.h"

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;

void Main(M& m, Vars& var) {
  auto sem = m.GetSem();
  struct {
    FieldNode<Scal> fnl;
    std::unique_ptr<Embed<M>> eb_;
    std::vector<std::pair<std::string, std::vector<Scal>>> data;
    FieldCell<Scal> fcu;
  } * ctx(sem);
  auto& eb_ = ctx->eb_;
  auto& fcu = ctx->fcu;
  if (sem("ctor")) {
    eb_.reset(new Embed<M>(m, var.Double["embed_gradlim"]));
  }
  if (sem.Nested("levelset")) {
    UEmbed<M>::InitLevelSet(ctx->fnl, m, var, false);
  }
  if (sem.Nested("init")) {
    eb_->Init(ctx->fnl);
  }
  if (sem.Nested()) {
    eb_->DumpPoly(var.Int["vtkbin"], var.Int["vtkmerge"]);
  }
  if (sem("dump")) {
    auto& eb = *eb_;
    fcu.Reinit(m);
    const Vect kA(-3., -3., 0.);
    // const Vect kA(0);
    const Vect kB(3., 3., 0.);
    const Vect kC(0.5, 0.5, 0.);
    auto func = [&](Vect x) { return (kA * x).dot(x) + kB.dot(x - kC); };
    auto funcg = [&](Vect x) { return kA * x * 2 + kB; };
    for (auto c : m.AllCells()) {
      fcu[c] = func(m.GetCenter(c));
    }
    std::vector<Scal> theta, grad, exact;

    if (m.IsRoot()) {
      auto c = m(m.FindNearestCell(Vect(0.4, 0.4, 0.)));
      const int imax = 100;
      for (int i = 0; i < imax; ++i) {
        const Scal t = 2 * M_PI * i / imax;
        const Vect nf(std::cos(t), std::sin(t), 0.);
        auto h = m.GetCellSize()[0];
        const Vect rf = c.center() + Vect(h * 0.2, h * 0.3, 0.);
        // const Scal g = GradDirichletQuad(rf, func(rf), nf, c, fcu, eb);
        const Scal g = GradDirichletQuadSecond(rf, func(rf), nf, c, fcu, eb);
        // const Scal g = GradDirichletLinear(rf, func(rf), nf, c, fcu, eb);
        // const Scal g = GradDirichletLinearFit(rf, func(rf), nf, c, fcu, eb);
        theta.push_back(t);
        grad.push_back(g);
        exact.push_back(funcg(rf).dot(nf));
      }
    }
    ctx->data.emplace_back(std::string("theta"), theta);
    ctx->data.emplace_back(std::string("grad"), grad);
    ctx->data.emplace_back(std::string("exact"), exact);
  }
  if (sem.Nested()) {
    DumpCsv(ctx->data, var.String["output_csv_omz"], m);
  }
  if (sem()) { // XXX empty stage
  }
}

int main(int argc, const char** argv) {
  MpiWrapper mpi(&argc, &argv);
  return RunMpiBasicFile<M>(mpi, Main);
}
