// Created by Petr Karnakov on 25.02.2020
// Copyright 2020 ETH Zurich

#include <iostream>

#include <util/posthook.h>

template <class M>
void PostHook(
    const Vars& var, const FieldCell<typename M::Vect>& fcvel, M& m,
    const Embed<M>& eb) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  auto sem = m.GetSem();
  struct {
    FieldCell<Scal> fc_theta;
    FieldCell<Scal> fc_theta_exact;
    FieldCell<Scal> fc_theta_error;
  } * ctx(sem);
  auto& fc_theta = ctx->fc_theta;
  auto& fc_theta_exact = ctx->fc_theta_exact;
  auto& fc_theta_error = ctx->fc_theta_error;
  if (sem()) {
    fc_theta.Reinit(eb, 0);
    fc_theta_error.Reinit(eb, 0);
    fc_theta_exact.Reinit(eb, 0);
    const Vect xc(0.5, 0.5, 0);
    const Scal r0 = 0.2;
    const Scal r1 = 0.4;

    for (auto c : eb.Cells()) {
      const Vect dx = m.GetCenter(c) - xc;
      const Vect dxuni = dx / dx.norm();
      fc_theta[c] = dxuni.cross_third(fcvel[c]);
      const Scal r = dx.norm();
      fc_theta_exact[c] = r * (sqr(r1 / r) - 1) / (sqr(r1 / r0) - 1) / r0;
      fc_theta_error[c] = fc_theta[c] - fc_theta_exact[c];
    }
    m.Dump(&fc_theta, "theta");
    m.Dump(&fc_theta_exact, "theta_exact");
    m.Dump(&fc_theta_error, "theta_error");
    m.Dump(&fcvel, 0, "vx");
    m.Dump(&fcvel, 1, "vy");
  }
  if (sem()) {
  }
  (void) var;
}

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using EB = Embed<M>;

template void PostHook(const Vars&, const FieldCell<Vect>&, M&, const EB&);
