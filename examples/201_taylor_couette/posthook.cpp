// Created by Petr Karnakov on 25.02.2020
// Copyright 2020 ETH Zurich

#include <fstream>

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
    Scal norm1, norm2, norminf, count;
  } * ctx(sem);
  auto& fc_theta = ctx->fc_theta;
  auto& fc_theta_exact = ctx->fc_theta_exact;
  auto& fc_theta_error = ctx->fc_theta_error;
  auto& norm1 = ctx->norm1;
  auto& norm2 = ctx->norm2;
  auto& norminf = ctx->norminf;
  auto& count = ctx->count;
  const Vect xc(0.5, 0.5, 0);
  const Scal r0 = 0.2;
  const Scal r1 = 0.4;
  if (sem()) {
    fc_theta.Reinit(eb, 0);
    fc_theta_error.Reinit(eb, 0);
    fc_theta_exact.Reinit(eb, 0);

    for (auto c : eb.Cells()) {
      Vect dx = m.GetCenter(c) - xc;
      dx[2] = 0;
      const Vect dxuni = dx / dx.norm();
      fc_theta[c] = dxuni.cross_third(fcvel[c]);
      const Scal r = dx.norm();
      fc_theta_exact[c] = r * (sqr(r1 / r) - 1) / (sqr(r1 / r0) - 1) / r0;
      fc_theta_error[c] = fc_theta[c] - fc_theta_exact[c];

      auto e = std::abs(fc_theta_error[c]);
      norm1 += e;
      norm2 += sqr(e);
      norminf = std::max(norminf, e);
      count += 1;
    }
    m.Dump(&fc_theta, "theta");
    m.Dump(&fc_theta_exact, "theta_exact");
    m.Dump(&fc_theta_error, "theta_error");
    m.Reduce(&norm1, "sum");
    m.Reduce(&norm2, "sum");
    m.Reduce(&norminf, "max");
    m.Reduce(&count, "sum");
  }
  if (sem()) {
    norm1 = norm1 / count;
    norm2 = std::sqrt(norm2 / count);
    if (m.IsRoot()) {
      std::ofstream fout("error");
      auto nx = m.GetGlobalSize()[0];
      auto cells_per_diameter = nx * r0 * 2;
      fout << cells_per_diameter << " ";
      fout << norm1 << " " << norm2 << " " << norminf << std::endl;
    }
  }
  (void) var;
}

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using EB = Embed<M>;

template void PostHook(const Vars&, const FieldCell<Vect>&, M&, const EB&);
