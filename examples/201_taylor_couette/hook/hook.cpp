// Created by Petr Karnakov on 25.02.2020
// Copyright 2020 ETH Zurich

#include <fstream>

#include <kernel/hydro.h>
#include <util/posthook.h>

// Local context, one for each block.
template <class M>
struct HookContext {
  using Scal = typename M::Scal;
  FieldCell<Scal> fc_theta;
  FieldCell<Scal> fc_theta_exact;
  FieldCell<Scal> fc_theta_error;
};

template <class M>
std::vector<generic::Vect<typename M::Scal, 3>> GetNorms(
    const std::vector<const FieldCell<typename M::Scal>*>& vfc, M& m,
    const Embed<M>& eb) {
  auto sem = m.GetSem(__func__);
  using Scal = typename M::Scal;
  using V = generic::Vect<Scal, 3>;
  struct {
    std::vector<V> vr;
    Scal vol;
  } * ctx(sem);
  auto& vr = ctx->vr;
  auto& vol = ctx->vol;
  if (sem("local")) {
    vr.resize(vfc.size(), V(0));
    for (size_t i = 0; i < vfc.size(); ++i) {
      auto& fc = *vfc[i];
      auto& r = vr[i];
      auto& norm1 = r[0];
      auto& norm2 = r[1];
      auto& norminf = r[2];
      for (auto c : eb.Cells()) {
        auto u = std::abs(fc[c]);
        norm1 += u * eb.GetVolume(c);
        norm2 += u * u * eb.GetVolume(c);
        norminf = std::max(norminf, u);
      }
      m.Reduce(&norm1, "sum");
      m.Reduce(&norm2, "sum");
      m.Reduce(&norminf, "max");
    }
    vol = 0;
    for (auto c : eb.Cells()) {
      (void)c;
      vol += eb.GetVolume(c);
    }
    m.Reduce(&vol, "sum");
  }
  if (sem("reduce")) {
    for (auto& r : vr) {
      auto& norm1 = r[0];
      auto& norm2 = r[1];
      norm1 = norm1 / vol;
      norm2 = std::sqrt(norm2 / vol);
    }
    return vr;
  }
  return {};
}

template <class M>
void InitHook(Hydro<M>* hydro) {
  auto& m = hydro->m;
  auto sem = m.GetSem();
  if (sem()) {
    // Create local context.
    hydro->hook_context_ =
        std::make_unique<Holder<HookContext<M>>>(new HookContext<M>());
  }
}

template <class M>
void PostHook(
    const Vars& var, const FieldCell<typename M::Vect>& fcvel, Hydro<M>* hydro,
    M& m, const Embed<M>& eb) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  auto sem = m.GetSem();
  struct {
    std::vector<generic::Vect<Scal, 3>> vnorms;
  } * ctx(sem);
  auto& t = *ctx;
  auto* hctx =
      dynamic_cast<Holder<HookContext<M>>*>(hydro->hook_context_.get())->Get();
  auto& ht = *hctx;
  const Vect xc(0.5, 0.5, 0);
  const Scal r0 = 0.2;
  const Scal r1 = 0.4;
  if (sem()) {
    auto* hctx =
        dynamic_cast<Holder<HookContext<M>>*>(hydro->hook_context_.get())->Get();
    ht.fc_theta.Reinit(eb, 0);
    ht.fc_theta_error.Reinit(eb, 0);
    ht.fc_theta_exact.Reinit(eb, 0);

    for (auto c : eb.Cells()) {
      Vect dx = m.GetCenter(c) - xc;
      dx[2] = 0;
      const Vect dxuni = dx / dx.norm();
      ht.fc_theta[c] = dxuni.cross_third(fcvel[c]);
      const Scal r = dx.norm();
      ht.fc_theta_exact[c] = r * (sqr(r1 / r) - 1) / (sqr(r1 / r0) - 1) / r0;
      ht.fc_theta_error[c] = ht.fc_theta[c] - ht.fc_theta_exact[c];
    }
    m.Dump(&ht.fc_theta, "theta");
    m.Dump(&ht.fc_theta_exact, "theta_exact");
    m.Dump(&ht.fc_theta_error, "theta_error");
  }
  if (sem.Nested()) {
    t.vnorms = GetNorms({&ht.fc_theta_error}, m, eb);
  }
  if (sem()) {
    if (m.IsRoot()) {
      std::ofstream fout("error");
      auto nx = m.GetGlobalSize()[0];
      auto cells_per_diameter = nx * r0 * 2;
      fout << cells_per_diameter << " ";
      const auto v = t.vnorms[0];
      fout << v[0] << " " << v[1] << " " << v[2] << std::endl;
    }
  }
}

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using EB = Embed<M>;

template void InitHook(Hydro<M>*);
template void PostHook(
    const Vars&, const FieldCell<Vect>&, Hydro<M>*, M&, const EB&);
