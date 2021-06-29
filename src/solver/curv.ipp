// Created by Petr Karnakov on 25.12.2019
// Copyright 2019 ETH Zurich

#include <array>
#include <exception>
#include <fstream>
#include <limits>
#include <memory>

#include "approx.h"
#include "debug/isnan.h"
#include "embed.h"
#include "geom/block.h"
#include "normal.h"
#include "parse/curv.h"
#include "reconst.h"
#include "trackerm.h"
#include "util/vof.h"
#include "vof.h"
#include "vofm.h"

#include "curv.h"

namespace curvature {

template <class M_>
struct Particles<M_>::Imp {
  Imp(M& m, const typename PartStrMeshM<M>::Par& par,
      const GRange<size_t>& layers)
      : partstrmeshm_(new PartStrMeshM<M>(m, par, layers)) {}

  template <class EB>
  void CalcCurvature(
      const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m,
      const EB& eb) {
    const auto& layers = plic.layers;
    auto sem = m.GetSem();

    if (sem.Nested("part")) {
      partstrmeshm_->Part(plic, eb);
    }
    if (sem("copy")) {
      fck.assert_size(layers);
      for (auto l : layers) {
        (*fck[l]) = *partstrmeshm_->GetCurv()[l];
      }
    }
  }

  std::unique_ptr<PartStrMeshM<M>> partstrmeshm_;
};

template <class M_>
Particles<M_>::Particles(
    M& m, const typename PartStrMeshM<M>::Par& par,
    const GRange<size_t>& layers)
    : imp(new Imp(m, par, layers)) {}

template <class EB_>
Particles<EB_>::~Particles() = default;

template <class M_>
void Particles<M_>::CalcCurvature(
    const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m, const M& eb) {
  imp->CalcCurvature(fck, plic, m, eb);
}

template <class M_>
void Particles<M_>::CalcCurvature(
    const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m,
    const Embed<M>& eb) {
  imp->CalcCurvature(fck, plic, m, eb);
}

template <class M_>
std::unique_ptr<PartStrMeshM<M_>> Particles<M_>::ReleaseParticles() {
  return std::move(imp->partstrmeshm_);
}

template <class M_>
const PartStrMeshM<M_>* Particles<M_>::GetParticles() const {
  return imp->partstrmeshm_.get();
}

template <class M_>
struct Heights<M_>::Imp {
  Imp() {}

  // fcu: volume fraction [a]
  // fcud2: volume fraction difference double (xpp-xmm, ypp-ymm, zpp-zmm) [i]
  static void CalcDiff2(
      const FieldCell<Scal>& fcu, FieldCell<Vect>& fcud2, const M& m) {
    fcud2.Reinit(m);
    for (auto c : m.CellsM()) {
      for (auto dr : m.dirs) {
        const auto d = m.direction(dr);
        const auto cmm = c - d - d;
        const auto cpp = c + d + d;
        fcud2[c][d] = fcu[cpp] - fcu[cmm];
      }
    }
  }
  // fcu: volume fraction [a]
  // fcud2: volume fraction difference double (xpp-xmm, ...) [a]
  // Output:
  // fcud4: volume fraction difference quad (xp4-xm4, ...) [i]
  static void CalcDiff4(
      const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcud2,
      FieldCell<Vect>& fcud4, const M& m) {
    fcud4.Reinit(m);
    for (auto c : m.CellsM()) {
      for (auto dr : m.dirs) {
        const auto d = m.direction(dr);
        const auto cmm = c - d - d;
        const auto cpp = c + d + d;
        const Scal um4 = fcu[c] - fcud2[cmm][d];
        const Scal up4 = fcud2[cpp][d] + fcu[c];
        fcud4[c][d] = up4 - um4;
      }
    }
  }
  // fcu: volume fraction [a]
  // fcud4: volume fraction difference double (xp4-xm4, ...) [a]
  // Output:
  // fcud6: volume fraction difference quad (xp6-xm6, ...) [i]
  static void CalcDiff6(
      const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcud4,
      FieldCell<Vect>& fcud6, const M& m) {
    fcud6.Reinit(m);
    for (auto c : m.CellsM()) {
      for (auto dr : m.dirs) {
        const auto d = m.direction(dr);
        const auto cmm = c - d - d;
        const auto cpp = c + d + d;
        const Scal um6 = fcu[cpp] - fcud4[cmm][d];
        const Scal up6 = fcud4[cpp][d] + fcu[cmm];
        fcud6[c][d] = up6 - um6;
      }
    }
  }

  template <class EB>
  void CalcCurvature(
      const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m, const EB&) {
    auto sem = m.GetSem();
    struct {
      FieldCell<Vect> fcud2; // difference of volume fractions with step 2
      FieldCell<Vect> fcud4; // difference of volume fractions with step 4
      FieldCell<Vect> fch; // height function
    } * ctx(sem);
    auto& t = *ctx;

    if (sem("diff2")) {
      Imp::CalcDiff2(*plic.vfcu[0], t.fcud2, m);
      m.Comm(&t.fcud2);
    }
    if (sem("diff4")) {
      Imp::CalcDiff4(*plic.vfcu[0], t.fcud2, t.fcud4, m);
      m.Comm(&t.fcud4);
    }
    if (sem("height")) {
      UNormal<M>::CalcHeight(
          m, *plic.vfcu[0], t.fcud2, t.fcud4, m.GetEdim(), t.fch);
      m.Comm(&t.fch);
    }
    if (sem("curvcomm")) {
      UNormal<M>::CalcCurvHeight(
          m, *plic.vfcu[0], t.fch, *plic.vfcn[0], m.GetEdim(), *fck[0]);
      m.Comm(fck[0]);
    }
  }
};

template <class M_>
Heights<M_>::Heights() = default;

template <class EB_>
Heights<EB_>::~Heights() = default;

template <class M_>
void Heights<M_>::CalcCurvature(
    const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m, const M& eb) {
  imp->CalcCurvature(fck, plic, m, eb);
}

template <class M_>
void Heights<M_>::CalcCurvature(
    const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m,
    const Embed<M>& eb) {
  imp->CalcCurvature(fck, plic, m, eb);
}

template <class M>
std::unique_ptr<Estimator<M>> MakeEstimator(
    const Vars& var, M& m, const GRange<size_t>& layers) {
  using Scal = typename M::Scal;
  const auto name = var.String("curvature", "particles");
  if (name == "particles") {
    const auto ps = ParsePar<PartStr<Scal>>()(m.GetCellSize()[0], var);
    const auto psm = ParsePar<PartStrMeshM<M>>()(ps, var);
    return std::make_unique<curvature::Particles<M>>(m, psm, layers);
  } else if (name == "heights") {
    return std::make_unique<curvature::Heights<M>>();
    //} else if (name == "hybrid") {
    //  return std::make_unique<curvature::Heights<M>>();
  }
  fassert(false, util::Format("Unknown curvature estimator '{}'", name));
}

} // namespace curvature
