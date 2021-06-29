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
#include "reconst.h"
#include "trackerm.h"
#include "util/vof.h"
#include "vof.h"
#include "vofm.h"

#include "curv.h"

template <class M_>
struct UCurv<M_>::Imp {
  static constexpr size_t dim = M::dim;

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
};

template <class M>
void UCurv<M>::CalcCurvHeight(
    const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcn, size_t edim,
    FieldCell<Scal>& fck, M& m) {
  auto sem = m.GetSem();
  struct {
    FieldCell<Vect> fcud2; // volume fraction difference double
    FieldCell<Vect> fcud4; // volume fraction difference quad
    FieldCell<Vect> fch; // height functions
  } * ctx(sem);
  auto& fcud2 = ctx->fcud2;
  auto& fcud4 = ctx->fcud4;
  auto& fch = ctx->fch;

  if (sem("diff2")) {
    Imp::CalcDiff2(fcu, fcud2, m);
    m.Comm(&fcud2);
  }
  if (sem("diff4")) {
    Imp::CalcDiff4(fcu, fcud2, fcud4, m);
    m.Comm(&fcud4);
  }
  if (sem("height")) {
    UNormal<M>::CalcHeight(m, fcu, fcud2, fcud4, edim, fch);
    m.Comm(&fch);
  }
  if (sem("curvcomm")) {
    UNormal<M>::CalcCurvHeight(m, fcu, fch, fcn, edim, fck);
    m.Comm(&fck);
  }
}

namespace curvature {

template <class M_>
struct Particles<M_>::Imp {
  Imp(const typename PartStrMeshM<M>::Par& par) : par_(par) {}

  template <class EB>
  void CalcCurvature(
      const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m,
      const EB& eb) {
    const auto& layers = plic.layers;
    auto sem = m.GetSem();

    if (sem("init")) {
      if (!psm_) {
        psm_.reset(new PartStrMeshM<M>(m, par_, layers));
      }
    }
    if (sem.Nested("part")) {
      psm_->Part(plic, eb);
    }
    if (sem("copy")) {
      fck.assert_size(layers);
      for (auto l : layers) {
        (*fck[l]) = *psm_->GetCurv()[l];
      }
    }
  }

  const typename PartStrMeshM<M>::Par& par_;
  std::unique_ptr<PartStrMeshM<M>> psm_;
};

template <class M_>
Particles<M_>::Particles(const typename PartStrMeshM<M>::Par& par)
    : imp(new Imp(par)) {}

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
  return std::move(imp->psm_);
}

template <class M_>
const PartStrMeshM<M_>* Particles<M_>::GetParticles() const {
  return imp->psm_.get();
}

} // namespace curvature
