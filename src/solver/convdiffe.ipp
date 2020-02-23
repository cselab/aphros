// Created by Petr Karnakov on 30.04.2019
// Copyright 2019 ETH Zurich

#include <cmath>
#include <sstream>
#include <stdexcept>

#include "approx.h"
#include "approx_eb.h"
#include "convdiffe.h"
#include "debug/isnan.h"

template <class EB_>
struct ConvDiffScalExp<EB_>::Imp {
  using Owner = ConvDiffScalExp<EB_>;
  using Vect = typename M::Vect;
  using UEB = UEmbed<M>;

  Imp(Owner* owner, const FieldCell<Scal>& fcu,
      const MapEmbed<BCond<Scal>>& mebc)
      : owner_(owner)
      , par(owner_->GetPar())
      , m(owner_->m)
      , eb(owner_->eb)
      , mebc_(mebc)
      , dtp_(-1.)
      , er_(0) {
    fcu_.time_curr = fcu;
    fcu_.time_prev = fcu;
    fcu_.iter_curr = fcu;
  }
  // Fields:
  void StartStep() {
    owner_->ClearIter();
    CHECKNAN(fcu_.time_curr, m.CN())

    fcu_.iter_curr = fcu_.time_curr;

    if (dtp_ == -1.) { // TODO: revise
      dtp_ = owner_->GetTimeStep();
    }

    er_ = 0.;
  }
  template <size_t dummy>
  FieldCell<Scal> CalcConvDiff(
      const FieldCell<Scal>& fcu, const FieldFace<Scal>& ffv) const {
    FieldCell<Scal> fclb(m, 0);

    {
      FieldFace<Scal> ffq; // convective fluxes
      const FieldCell<Vect> fcg = Gradient(UEB::Interpolate(fcu, mebc_, m), m);
      ffq = UEB::InterpolateUpwind(fcu, mebc_, par.sc, fcg, ffv, m);
      for (auto f : m.Faces()) {
        ffq[f] *= (*owner_->ffv_)[f];
      }
      for (auto c : m.Cells()) {
        Scal sum = 0.; // sum
        for (auto q : m.Nci(c)) {
          const IdxFace f = m.GetFace(c, q);
          sum += ffq[f] * m.GetOutwardFactor(c, q);
        }
        fclb[c] += sum / m.GetVolume(c) * (*owner_->fcr_)[c];
      }
    }

    if (owner_->ffd_) {
      FieldFace<Scal> ffq; // diffusive fluxes
      ffq = UEB::Gradient(fcu, mebc_, m);
      for (auto f : m.Faces()) {
        ffq[f] *= (-(*owner_->ffd_)[f]) * m.GetArea(f);
      }
      for (auto c : m.Cells()) {
        Scal sum = 0.; // sum
        for (auto q : m.Nci(c)) {
          const IdxFace f = m.GetFace(c, q);
          sum += ffq[f] * m.GetOutwardFactor(c, q);
        }
        fclb[c] += sum / m.GetVolume(c);
      }
    }
    return fclb;
  }
  template <size_t dummy>
  FieldCell<Scal> CalcConvDiff(
      const FieldCell<Scal>& fcu, const FieldEmbed<Scal>& fev) const {
    const FieldCell<Vect> fcg = Gradient(UEB::Interpolate(fcu, mebc_, m), m);

    FieldCell<Scal> fclb(m, 0);
    {
      // convective fluxes
      FieldEmbed<Scal> feq =
          UEB::InterpolateUpwind(fcu, mebc_, par.sc, fcg, fev, eb);
      for (auto f : eb.Faces()) {
        feq[f] *= fev[f];
      }
      for (auto c : eb.CFaces()) {
        feq[c] *= fev[c];
      }
      for (auto c : eb.Cells()) {
        Scal sum = feq[c];
        for (auto q : eb.Nci(c)) {
          const IdxFace f = m.GetFace(c, q);
          sum += feq[f] * m.GetOutwardFactor(c, q);
        }
        fclb[c] += sum * (*owner_->fcr_)[c];
      }
    }

    if (owner_->ffd_) {
      // diffusive fluxes
      FieldEmbed<Scal> feq = UEB::Gradient(fcu, mebc_, eb);
      for (auto f : eb.Faces()) {
        feq[f] *= -(*owner_->ffd_)[f] * eb.GetArea(f);
      }
      for (auto c : eb.CFaces()) {
        feq[c] *= -(*owner_->ffd_)[c] * eb.GetArea(c);
      }
      for (auto c : eb.Cells()) {
        Scal sum = feq[c];
        for (auto q : eb.Nci(c)) {
          IdxFace f = m.GetFace(c, q);
          sum += feq[f] * m.GetOutwardFactor(c, q);
        }
        fclb[c] += sum;
      }
    }

    fclb = UEB::RedistributeCutCells(fclb, eb);
    for (IdxCell c : eb.Cells()) {
      fclb[c] /= eb.GetVolume(c);
    }
    return fclb;
  }
  // Assembles linear system
  // fcu: field from previous iteration [a]
  // ffv: volume flux
  // Output:
  // diagonal linear system: fcla * (fcup - fcu) + fclb = 0
  // where fcup at new iteration
  void Assemble(
      const FieldCell<Scal>& fcu, const FieldFaceb<Scal>& ffv,
      FieldCell<Scal>& fcla, FieldCell<Scal>& fclb) const {
    auto sem = m.GetSem("assemble");

    if (sem("assemble")) {
      const FieldCell<Vect> fcg = Gradient(UEB::Interpolate(fcu, mebc_, m), m);

      fcla.Reinit(m);
      fclb = CalcConvDiff<0>(fcu, ffv);

      // time derivative coeffs
      const Scal dt = owner_->GetTimeStep();
      const std::vector<Scal> ac =
          GetGradCoeffs(0., {-(dt + dtp_), -dt, 0.}, par.second ? 0 : 1);

      for (IdxCell c : eb.Cells()) {
        // time derivative
        Scal r = (*owner_->fcr_)[c];
        fcla[c] = ac[2] * r;
        fclb[c] += (ac[0] * fcu_.time_prev[c] + ac[1] * fcu_.time_curr[c]) * r;
        // source
        fclb[c] += -(*owner_->fcs_)[c];
        // delta form
        fclb[c] += fcla[c] * fcu[c];
        // under-relaxation
        fcla[c] /= par.relax;
      }
    }
  }
  // Assembles linear system
  // fcu: field from previous iteration [a]
  // ffv: volume flux
  // Output:
  // fcl: linear system
  void Assemble(const FieldCell<Scal>& fcu, const FieldFaceb<Scal>& ffv) {
    Assemble(fcu, ffv, fcla_, fclb_);
  }
  // Solves linear system.
  // fcla * u + fclb = 0: system to solve
  // Output:
  // fcu: result
  void Solve(
      const FieldCell<Scal>& fcla, const FieldCell<Scal>& fclb,
      FieldCell<Scal>& fcu) {
    fcu.Reinit(m, 0);
    for (auto c : eb.Cells()) {
      fcu[c] = -fclb[c] / fcla[c];
    }
  }
  void MakeIteration() {
    auto sem = m.GetSem("convdiff-iter");

    auto& prev = fcu_.iter_prev;
    auto& curr = fcu_.iter_curr;

    if (sem("init")) {
      prev.swap(curr);
    }
    if (sem.Nested("assemble")) {
      Assemble(prev, *owner_->ffv_, fcla_, fclb_);
    }
    if (sem("solve")) {
      Solve(fcla_, fclb_, curr); // solve for correction, store in curr
    }
    if (sem("apply")) {
      CalcError();

      // apply, store result in curr
      for (auto c : eb.Cells()) {
        curr[c] += prev[c];
      }
      m.Comm(&curr);
      owner_->IncIter();
    }
  }
  void FinishStep() {
    fcu_.time_prev.swap(fcu_.time_curr);
    fcu_.time_curr = fcu_.iter_curr;
    CHECKNAN(fcu_.time_curr, m.CN())
    owner_->IncTime();
    dtp_ = owner_->GetTimeStep();
  }
  void CalcError() {
    // max difference between iter_curr and iter_prev
    auto& prev = fcu_.iter_prev;
    auto& curr = fcu_.iter_curr;
    Scal a = 0;
    for (auto c : eb.Cells()) {
      a = std::max<Scal>(a, std::abs(curr[c] - prev[c]));
    }
    er_ = a;
  }
  double GetError() const {
    return er_;
  }
  // Apply correction to field and comm
  // uc: correction [i]
  // Output:
  // u(l) += uc [a]
  void CorrectField(Step l, const FieldCell<Scal>& uc) {
    auto sem = m.GetSem("corr");
    if (sem("apply")) {
      auto& u = fcu_.Get(l);
      for (auto c : eb.Cells()) {
        u[c] += uc[c];
      }
      CalcError();
      m.Comm(&u);
    }
  }
  FieldCell<Scal> GetDiag() const {
    return fcla_;
  }
  FieldCell<Scal> GetConst() const {
    return fclb_;
  }

  Owner* owner_;
  Par par;
  M& m; // mesh
  const EB& eb; // embed mesh

  StepData<FieldCell<Scal>> fcu_; // field
  const MapEmbed<BCond<Scal>>& mebc_; // boundary conditions

  size_t bc_ = 0; // boundary conditions, 0: value, 1: gradient
  Scal bcu_ = 0.; // value or grad.dot.outer_normal
  // diagonal linear system: fcla * (fcup - fcu) + fclb = 0
  FieldCell<Scal> fcla_;
  FieldCell<Scal> fclb_;

  Scal dtp_; // dt prev
  Scal er_; // error
};

template <class EB_>
ConvDiffScalExp<EB_>::ConvDiffScalExp(
    M& m, const EB& eb, const FieldCell<Scal>& fcu,
    const MapEmbed<BCond<Scal>>& mebc, const FieldCell<Scal>* fcr,
    const FieldFaceb<Scal>* ffd, const FieldCell<Scal>* fcs,
    const FieldFaceb<Scal>* ffv, double t, double dt, Par par)
    : Base(t, dt, m, eb, par, fcr, ffd, fcs, ffv)
    , imp(new Imp(this, fcu, mebc)) {}

template <class EB_>
ConvDiffScalExp<EB_>::~ConvDiffScalExp() = default;

template <class EB_>
void ConvDiffScalExp<EB_>::Assemble(
    const FieldCell<Scal>& fcu, const FieldFaceb<Scal>& ffv) {
  imp->Assemble(fcu, ffv);
}

template <class EB_>
void ConvDiffScalExp<EB_>::CorrectField(Step l, const FieldCell<Scal>& uc) {
  imp->CorrectField(l, uc);
}

template <class EB_>
auto ConvDiffScalExp<EB_>::GetDiag() const -> FieldCell<Scal> {
  return imp->GetDiag();
}

template <class EB_>
auto ConvDiffScalExp<EB_>::GetConst() const -> FieldCell<Scal> {
  return imp->GetConst();
}

template <class EB_>
void ConvDiffScalExp<EB_>::StartStep() {
  imp->StartStep();
}

template <class EB_>
void ConvDiffScalExp<EB_>::MakeIteration() {
  imp->MakeIteration();
}

template <class EB_>
void ConvDiffScalExp<EB_>::FinishStep() {
  imp->FinishStep();
}

template <class EB_>
double ConvDiffScalExp<EB_>::GetError() const {
  return imp->GetError();
}

template <class EB_>
auto ConvDiffScalExp<EB_>::GetField(Step l) const -> const FieldCell<Scal>& {
  return imp->fcu_.Get(l);
}
