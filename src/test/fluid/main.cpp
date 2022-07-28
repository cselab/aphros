// Created by Petr Karnakov on 03.02.2020
// Copyright 2020 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <random>
#include <sstream>
#include <typeinfo>

#include "distr/distrbasic.h"
#include "dump/dump.h"
#include "dump/hdf.h"
#include "geom/mesh.h"
#include "solver/approx.h"
#include "solver/approx_eb.h"
#include "solver/cond.h"
#include "solver/embed.h"
#include "solver/solver.h"
#include "util/convdiff.h"
#include "util/fluid.h"

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using EB = Embed<M>;
using UEB = UEmbed<M>;

template <class MEB>
void Test(M& m, MEB& eb, std::string name) {
  auto sem = m.GetSem();
  struct {
    FieldCell<Vect> fc_force;
    FieldCell<Vect> fc_vel;
  } * ctx(sem);
  auto& t = *ctx;
  if (sem()) {
    std::cout << name << " " << typeid(MEB).name() << std::endl;
    const Vect v0(0);
    const Vect vx(1., 2., 3.);
    const Vect vy = vx * 2;
    const Vect vz = vx * 3;

    t.fc_vel.Reinit(m, Vect(0));
    FieldFace<Scal> ff_mu(eb);
    for (auto c : eb.AllCells()) {
      auto x = m.GetCenter(c);
      t.fc_vel[c] = v0 + vx * x[0] + vy * x[1] + vz * x[2];
    }
    for (auto f : eb.Faces()) {
      auto x = eb.GetFaceCenter(f);
      x = x * x;
      ff_mu[f] = x[0] + x[1] * 10 + x[2] * 100;
    }

    t.fc_force.Reinit(m, Vect(0));
    UFluid<M>::AppendExplViscous(t.fc_force, t.fc_vel, {}, ff_mu, eb);
    m.Dump(&t.fc_force, 0, "fx");
    m.Dump(&t.fc_force, 1, "fy");
    m.Dump(&t.fc_force, 2, "fz");
    m.DumpCommit();
  }
  if (sem()) {
  }
}

void Main(M& m, Vars& var) {
  auto sem = m.GetSem();
  struct {
    FieldNode<Scal> fnl;
    std::unique_ptr<Embed<M>> eb_;
  } * ctx(sem);
  auto& eb_ = ctx->eb_;
  if (sem("ctor")) {
    eb_.reset(new Embed<M>(m, 0));
  }
  if (sem.Nested("levelset")) {
    UEmbed<M>::InitLevelSet(ctx->fnl, m, var, false);
  }
  if (sem.Nested("init")) {
    eb_->Init(ctx->fnl);
  }
  if (sem.Nested()) {
    eb_->DumpPoly(true, true);
  }
  if (sem.Nested()) {
    Test(m, *eb_, "mesh");
  }
}

int main(int argc, const char** argv) {
  MpiWrapper mpi(&argc, &argv);
  return RunMpiBasicFile<M>(mpi, Main);
}
