// Created by Petr Karnakov on 05.07.2020
// Copyright 2020 ETH Zurich

#undef NDEBUG
#include <mpi.h>
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
#include "solver/embed.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;

void Main(M& m, Vars& var) {
  auto sem = m.GetSem();
  struct {
    FieldNode<Scal> fnl;
    std::unique_ptr<Embed<M>> eb_;
    std::vector<std::pair<std::string, std::vector<Scal>>> data;
  } * ctx(sem);
  auto& eb_ = ctx->eb_;
  if (sem("ctor")) {
    eb_.reset(new Embed<M>(m, var.Double["embed_gradlim"]));
    ctx->fnl = UEmbed<M>::InitEmbed(m, var, false);
  }
  if (sem.Nested("init")) {
    eb_->Init(ctx->fnl);
  }
  if (sem("dump")) {
    auto& eb = *eb_;
    std::vector<Scal> vx, vy;
    for (auto c : eb.CFaces()) {
      vx.push_back(eb.GetFaceCenter(c)[0]);
      vy.push_back(eb.GetFaceCenter(c)[1]);
    }
    ctx->data.emplace_back(std::string("vx"), vx);
    ctx->data.emplace_back(std::string("vy"), vy);
  }
  if (sem.Nested()) {
    DumpCsv(ctx->data, var.String["output_csv_omz"], m);
  }
  if (sem()) { // XXX empty stage
  }
}

int main(int argc, const char** argv) {
  return RunMpiBasicFile<M>(argc, argv, Main);
}
