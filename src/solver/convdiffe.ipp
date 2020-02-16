// Created by Petr Karnakov on 30.04.2019
// Copyright 2019 ETH Zurich

#include <cmath>
#include <sstream>
#include <stdexcept>

#include "approx.h"
#include "convdiffe.h"
#include "debug/isnan.h"

template <class EB_>
struct ConvDiffScalExp<EB_>::Imp {
  using Owner = ConvDiffScalExp<EB_>;
  using Vect = typename M::Vect;

  Imp(Owner* owner, const EB& eb0, const FieldCell<Scal>& fcu,
      const MapCondFace& mfc)
      : owner_(owner)
      , par(owner_->GetPar())
      , m(owner_->m)
      , eb(eb0)
      , mfc_(mfc)
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
  // Assembles linear system
  // fcu: field from previous iteration [a]
  // ffv: volume flux
  // Output:
  // diagonal linear system: fcla * (fcup - fcu) + fclb = 0
  // where fcup at new iteration
  void Assemble(
      const FieldCell<Scal>& fcu, const FieldFaceb<Scal>& ffv,
      FieldCell<Scal>& fcla, FieldCell<Scal>& fclb) {
    auto sem = m.GetSem("assemble");

    if (sem("assemble")) {
      FieldCell<Vect> fcg = Gradient(Interpolate(fcu, mfc_, m), m);

      fcla.Reinit(m);
      fclb.Reinit(m, 0.);

      FieldFaceb<Scal> ffq; // flux tmp

      // Calc convective fluxes
      Interpolate(fcu, fcg, mfc_, ffv, m, par.sc, par.th, ffq);
      for (auto f : m.Faces()) {
        ffq[f] *= (*owner_->ffv_)[f];
      }
      // Append
      for (IdxCell c : m.Cells()) {
        Scal s = 0.; // sum
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetFace(c, q);
          s += ffq[f] * m.GetOutwardFactor(c, q);
        }
        fclb[c] += s / m.GetVolume(c) * (*owner_->fcr_)[c];
      }

      if (owner_->ffd_) {
        // Calc diffusive fluxes
        Gradient(fcu, mfc_, m, ffq);
        for (auto f : m.Faces()) {
          ffq[f] *= (-(*owner_->ffd_)[f]) * m.GetArea(f);
        }
        // Append
        for (IdxCell c : m.Cells()) {
          Scal s = 0.; // sum
          for (auto q : m.Nci(c)) {
            IdxFace f = m.GetFace(c, q);
            s += ffq[f] * m.GetOutwardFactor(c, q);
          }
          fclb[c] += s / m.GetVolume(c);
        }
      }

      // time derivative coeffs
      Scal dt = owner_->GetTimeStep();
      std::vector<Scal> ac =
          GetGradCoeffs(0., {-(dt + dtp_), -dt, 0.}, par.second ? 0 : 1);

      for (IdxCell c : m.Cells()) {
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
    fcu.Reinit(m);
    for (auto c : m.Cells()) {
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
      for (auto c : m.Cells()) {
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
    for (auto c : m.Cells()) {
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
      for (auto c : m.Cells()) {
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
  const MapCondFace& mfc_; // face cond

  // diagonal linear system: fcla * (fcup - fcu) + fclb = 0
  FieldCell<Scal> fcla_;
  FieldCell<Scal> fclb_;

  Scal dtp_; // dt prev
  Scal er_; // error
};

template <class EB_>
ConvDiffScalExp<EB_>::ConvDiffScalExp(
    M& m, const EB& eb, const FieldCell<Scal>& fcu, const MapCondFace& mfc,
    const FieldCell<Scal>* fcr, const FieldFaceb<Scal>* ffd,
    const FieldCell<Scal>* fcs, const FieldFaceb<Scal>* ffv, double t,
    double dt, Par par)
    : ConvDiffScal<M>(t, dt, m, par, fcr, ffd, fcs, ffv)
    , imp(new Imp(this, eb, fcu, mfc)) {}

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
