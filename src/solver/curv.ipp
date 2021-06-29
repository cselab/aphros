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
#include "util/height.h"
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
  // fcdu2: volume fraction difference double (xpp-xmm, ypp-ymm, zpp-zmm) [i]
  static void CalcDiff2(
      const FieldCell<Scal>& fcu, FieldCell<Vect>& fcdu2, const M& m) {
    fcdu2.Reinit(m);
    for (auto c : m.CellsM()) {
      for (auto dr : m.dirs) {
        const auto d = m.direction(dr);
        const auto cmm = c - d - d;
        const auto cpp = c + d + d;
        fcdu2[c][d] = fcu[cpp] - fcu[cmm];
      }
    }
  }
  // fcu: volume fraction [a]
  // fcdu2: volume fraction difference double (xpp-xmm, ...) [a]
  // Output:
  // fcdu4: volume fraction difference quad (xp4-xm4, ...) [i]
  static void CalcDiff4(
      const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcdu2,
      FieldCell<Vect>& fcdu4, const M& m) {
    fcdu4.Reinit(m);
    for (auto c : m.CellsM()) {
      for (auto dr : m.dirs) {
        const auto d = m.direction(dr);
        const auto cmm = c - d - d;
        const auto cpp = c + d + d;
        const Scal um4 = fcu[c] - fcdu2[cmm][d];
        const Scal up4 = fcdu2[cpp][d] + fcu[c];
        fcdu4[c][d] = up4 - um4;
      }
    }
  }
  // fcu: volume fraction [a]
  // fcdu4: volume fraction difference double (xp4-xm4, ...) [a]
  // Output:
  // fcdu6: volume fraction difference quad (xp6-xm6, ...) [i]
  static void CalcDiff6(
      const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcdu4,
      FieldCell<Vect>& fcdu6, const M& m) {
    fcdu6.Reinit(m);
    for (auto c : m.CellsM()) {
      for (auto dr : m.dirs) {
        const auto d = m.direction(dr);
        const auto cmm = c - d - d;
        const auto cpp = c + d + d;
        const Scal um6 = fcu[cpp] - fcdu4[cmm][d];
        const Scal up6 = fcdu4[cpp][d] + fcu[cmm];
        fcdu6[c][d] = up6 - um6;
      }
    }
  }

  // Computes heights.
  // S: the number of stages for stencils, column size is [-S*2,S*2]
  // fcu: volume fraction [a]
  // fcdu2: difference of volume fractions with step 2 [a]
  // fcdu4: difference of volume fractions with step 4 [a]
  // Output:
  // fch: fch[c][d] is absolute position of the interface
  // from a column in direction d starting from an interfacial cell c
  // otherwise, NaN
  template <size_t S>
  static void CalcHeight(
      M& m, const FieldCell<Scal>& fcu,
      const std::array<const FieldCell<Vect>*, S> vfcud, FieldCell<Vect>& fch) {
    auto I = [](Scal a) { return a > 0 && a < 1; }; // interface

    fch.Reinit(m, GetNan<Vect>());

    for (auto c : m.CellsM()) {
      if (!I(fcu[c])) {
        continue;
      }
      for (size_t dr = 0; dr < m.GetEdim(); ++dr) {
        const auto d = m.direction(dr);
        const auto cm = c - d;
        const auto cmm = c - d - d;
        const auto cp = c + d;
        const auto cpp = c + d + d;

        const size_t si = (S + 1) * 4 + 1;
        const size_t sih = (S + 1) * 2;

        std::array<Scal, si> uu;

        uu[sih] = fcu[c];
        uu[sih - 1] = fcu[cm];
        uu[sih + 1] = fcu[cp];
        uu[sih - 2] = fcu[cmm];
        uu[sih + 2] = fcu[cpp];

        for (size_t s = 0; s < S; ++s) {
          const size_t q = (s + 1) * 2;
          const size_t ia = q + 1;
          const size_t ib = q + 2;
          const FieldCell<Vect>& fcdu = *vfcud[s];

          uu[sih - ia] = uu[sih - ia + 2 * q] - fcdu[cm][d];
          uu[sih + ia] = uu[sih + ia - 2 * q] + fcdu[cp][d];
          uu[sih - ib] = uu[sih - ib + 2 * q] - fcdu[cmm][d];
          uu[sih + ib] = uu[sih + ib - 2 * q] + fcdu[cpp][d];
        }

        // |cm6|cm5|cm4|cm3|cmm| cm| c |cp |cpp|cp3|cp4|cp5|cp6|
        // |   |   |   |   | * |   | c |   | * |   |   |   |   |
        // |   |   | * |   |   |   | c |   |   |   | * |   |   |

        const Scal s = UHeight<Scal>::Good(uu);
        fch[c][d] = m.GetCenter(c)[d] + s * m.GetCellSize()[d];
      }
    }
  }

  // Computes curvature from height functions.
  // fcu: volume fraction
  // fcn: normal, antigradient of fcu
  // Output: modified in cells with fci=1, resized to m
  // fck: curvature [i]
  static void CalcCurvHeight(
      M& m, const FieldCell<Scal>& fcu, const FieldCell<Vect>& fch,
      const FieldCell<Vect>& fcn, FieldCell<Scal>& fck) {
    using MIdx = typename M::MIdx;
    using Dir = typename M::Dir;
    auto& bc = m.GetIndexCells();

    auto I = [](Scal a) { return a > 0 && a < 1; }; // interface

    fck.Reinit(m, GetNan<Scal>());

    for (auto c : m.Cells()) {
      if (!I(fcu[c])) {
        continue;
      }

      size_t di = fcn[c].abs().argmax(); // best direction index
      Dir dn(di); // best direction
      // directions of plane tangents ([d]irection [t]angents)
      Dir dtx((size_t(dn) + 1) % M::dim);
      Dir dty((size_t(dn) + 2) % M::dim);

      MIdx w = bc.GetMIdx(c);

      // offset in normal direction
      MIdx on = MIdx(dn);
      // offset in dtx,dty
      MIdx otx = MIdx(dtx);
      MIdx oty = MIdx(dty);
      // mesh step
      const Scal lx = m.GetCellSize()[size_t(dtx)];
      const Scal ly = m.GetCellSize()[size_t(dty)];

      // Evaluates height function from nearest interface
      // o: offset from w
      auto hh = [I, &w, &fch, &di, &bc, &on, &fcu](MIdx o) -> Scal {
        const int si = 5; // stencil size
        const int sih = si / 2;
        int i = sih; // closest interface to center
        while (i < si) {
          auto cn = bc.GetIdx(w + o + on * (i - sih));
          if (I(fcu[cn])) {
            return fch[cn][di];
          }
          if (i > sih) {
            i = si - i - 1;
          } else {
            i = si - i;
          }
        }
        return GetNan<Scal>();
      };

      // height function
      const Scal hcc = hh(MIdx(0));
      const Scal hmc = hh(-otx);
      const Scal hpc = hh(otx);
      const Scal hcm = hh(-oty);
      const Scal hcp = hh(oty);
      // corners: hxy
      const Scal hmm = hh(-otx - oty);
      const Scal hmp = hh(-otx + oty);
      const Scal hpm = hh(otx - oty);
      const Scal hpp = hh(otx + oty);

      // first derivative (slope)
      Scal hx = (hpc - hmc) / (2. * lx); // centered
      Scal hy = (hcp - hcm) / (2. * ly);
      // second derivative
      const Scal fl = 0.2; // filter factor (Basilisk: fl=0.2)
      Scal hxx = ((hpm - 2. * hcm + hmm) * fl + (hpc - 2. * hcc + hmc) +
                  (hpp - 2. * hcp + hmp) * fl) /
                 ((1 + 2 * fl) * lx * lx);
      Scal hyy = ((hmp - 2. * hmc + hmm) * fl + (hcp - 2. * hcc + hcm) +
                  (hpp - 2. * hpc + hpm) * fl) /
                 ((1 + 2 * fl) * ly * ly);
      Scal hxy = ((hpp - hmp) - (hpm - hmm)) / (4. * lx * ly);
      // curvature
      Scal k =
          (2. * hx * hy * hxy - (sqr(hy) + 1.) * hxx - (sqr(hx) + 1.) * hyy) /
          std::pow(sqr(hx) + sqr(hy) + 1., 3. / 2.);

      if (fcn[c][di] < 0) {
        k *= -1;
      }

      // curvature
      fck[c] = k;
    }
  }

  template <class EB>
  void CalcCurvature(
      const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m, const EB&) {
    auto sem = m.GetSem();
    struct {
      FieldCell<Vect> fcdu2; // difference of volume fractions with step 2
      FieldCell<Vect> fcdu4; // difference of volume fractions with step 4
      FieldCell<Vect> fch; // height function
    } * ctx(sem);
    auto& t = *ctx;

    if (sem("diff2")) {
      CalcDiff2(*plic.vfcu[0], t.fcdu2, m);
      m.Comm(&t.fcdu2);
    }
    if (sem("diff4")) {
      CalcDiff4(*plic.vfcu[0], t.fcdu2, t.fcdu4, m);
      m.Comm(&t.fcdu4);
    }
    if (sem("height")) {
      std::array<const FieldCell<Vect>*, 2> vfcdu = {&t.fcdu2, &t.fcdu4};
      CalcHeight(m, *plic.vfcu[0], vfcdu, t.fch);
      m.Comm(&t.fch);
    }
    if (sem("curvcomm")) {
      CalcCurvHeight(m, *plic.vfcu[0], t.fch, *plic.vfcn[0], *fck[0]);
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
