// Created by Petr Karnakov on 25.12.2019
// Copyright 2019 ETH Zurich

#include <array>
#include <exception>
#include <fstream>
#include <limits>
#include <memory>

#include "approx.h"
#include "debug/isnan.h"
#include "geom/block.h"
#include "normal.h"
#include "reconst.h"
#include "trackerm.h"
#include "tvd.h"
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
    for (auto c : m.Cells()) {
      for (size_t d = 0; d < dim; ++d) {
        auto cm = m.GetCell(c, 2 * d);
        auto cmm = m.GetCell(cm, 2 * d);
        auto cp = m.GetCell(c, 2 * d + 1);
        auto cpp = m.GetCell(cp, 2 * d + 1);
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
    for (auto c : m.Cells()) {
      for (size_t d = 0; d < dim; ++d) {
        auto cm = m.GetCell(c, 2 * d);
        auto cmm = m.GetCell(cm, 2 * d);
        auto cp = m.GetCell(c, 2 * d + 1);
        auto cpp = m.GetCell(cp, 2 * d + 1);
        Scal um4 = fcu[c] - fcud2[cmm][d];
        Scal up4 = fcud2[cpp][d] + fcu[c];
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
    for (auto c : m.Cells()) {
      for (size_t d = 0; d < dim; ++d) {
        auto cm = m.GetCell(c, 2 * d);
        auto cmm = m.GetCell(cm, 2 * d);
        auto cp = m.GetCell(c, 2 * d + 1);
        auto cpp = m.GetCell(cp, 2 * d + 1);
        Scal um6 = fcu[cpp] - fcud4[cmm][d];
        Scal up6 = fcud4[cpp][d] + fcu[cmm];
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

template <class M>
std::unique_ptr<PartStrMeshM<M>> UCurv<M>::CalcCurvPart(
    const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fca,
    const Multi<const FieldCell<Vect>*>& fcn,
    const Multi<const FieldCell<bool>*>& fci,
    const Multi<const FieldCell<Scal>*>& fccl,
    const typename PartStrMeshM<M>::Par& par,
    const Multi<FieldCell<Scal>*>& fck, M& m) {
  using PSM = PartStrMeshM<M>;

  auto sem = m.GetSem();
  struct {
    std::unique_ptr<PSM> psm;
  } * ctx(sem);
  auto& psm = ctx->psm;

  if (sem("init")) {
    psm.reset(new PSM(m, par, layers));
  }
  if (sem.Nested("part")) {
    psm->Part(fca, fcn, fci, fccl);
  }
  if (sem("copy")) {
    fck.assert_size(layers);
    for (auto l : layers) {
      (*fck[l]) = *psm->GetCurv()[l];
    }
    return std::move(psm);
  }
  return nullptr;
}

template <class M>
std::unique_ptr<PartStrMeshM<M>> UCurv<M>::CalcCurvPart(
    const GRange<size_t>& layers, const AdvectionSolver<M>* asbase,
    const typename PartStrMeshM<M>::Par& par,
    const Multi<FieldCell<Scal>*>& fck, M& m) {
  if (auto as = dynamic_cast<const Vof<M>*>(asbase)) {
    return CalcCurvPart(
        layers, &as->GetAlpha(), &as->GetNormal(), &as->GetMask(), nullptr, par,
        fck, m);
  } else if (auto as = dynamic_cast<const Vofm<M>*>(asbase)) {
    return CalcCurvPart(
        layers, as->GetAlpha(), as->GetNormal(), as->GetMask(), as->GetColor(),
        par, fck, m);
  } else {
    for (auto l : layers) {
      fck[l]->Reinit(m, GetNan<Scal>());
    }
    return nullptr;
  }
}
