// Created by Petr Karnakov on 30.07.2018
// Copyright 2018 ETH Zurich

#include <cmath>
#include <sstream>
#include <stdexcept>

#include "approx.h"
#include "approx_eb.h"
#include "convdiffi.h"
#include "debug/isnan.h"
#include "linear/linear.h"
#include "util/convdiff.h"

template <class EB_>
struct ConvDiffScalImp<EB_>::Imp {
  using Owner = ConvDiffScalImp<EB_>;
  using Vect = typename M::Vect;
  using ExprFace = typename M::ExprFace;
  using Expr = typename M::Expr;
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
  Scal Eval(const Expr& e, IdxCell c, const FieldCell<Scal>& fcu) {
    Scal r = e[Expr::dim - 1];
    r += fcu[c] * e[0];
    for (auto q : eb.Nci(c)) {
      r += fcu[eb.GetCell(c, q)] * e[1 + q];
    }
    return r;
  }
  // Assembles linear system
  // fcu: field from previous iteration [a]
  // ffv: volume flux
  // Output:
  // fcl: linear system
  void Assemble(
      const FieldCell<Scal>& fcu, const FieldFaceb<Scal>& ffv,
      FieldCell<Expr>& fcl) {
    FieldFaceb<ExprFace> ffq(eb, ExprFace(0));
    // convective fluxes
    {
      const FieldCell<Vect> fcg =
          UEB::Gradient(UEB::Interpolate(fcu, mebc_, eb), eb);
      const FieldFaceb<ExprFace> ffu = UEB::InterpolateUpwindImplicit(
          fcu, mebc_, par.sc, par.df, fcg, ffv, eb);
      eb.LoopFaces([&](auto cf) { //
        ffq[cf] += ffu[cf] * ffv[cf];
      });
    }
    // diffusive fluxes
    if (owner_->ffd_) {
      const FieldFaceb<ExprFace> ffg = UEB::GradientImplicit(fcu, mebc_, eb);
      eb.LoopFaces([&](auto cf) { //
        ffq[cf] -= ffg[cf] * (*owner_->ffd_)[cf] * eb.GetArea(cf);
      });
    }

    const Scal dt = owner_->GetTimeStep();
    std::vector<Scal> time_coeff =
        GetGradCoeffs(0., {-(dt + dtp_), -dt, 0.}, par.second ? 0 : 1);

    fcl.Reinit(eb, Expr::GetUnit(0)); // initialize as diagonal system
    for (auto c : eb.Cells()) {
      Expr sum(0);
      eb.LoopNci(c, [&](auto q) {
        const auto cf = eb.GetFace(c, q);
        const ExprFace v = ffq[cf] * eb.GetOutwardFactor(c, q);
        eb.AppendExpr(sum, v, q);
      });

      Expr td(0); // time derivative
      td[0] = time_coeff[2];
      td[Expr::dim - 1] =
          time_coeff[0] * fcu_.time_prev[c] + time_coeff[1] * fcu_.time_curr[c];

      Expr& e = fcl[c];
      e = (td + sum / eb.GetVolume(c)) * (*owner_->fcr_)[c];
      e[Expr::dim - 1] -= (*owner_->fcs_)[c];

      // Convert to delta-form
      e[Expr::dim - 1] = Eval(e, c, fcu);
      // Apply under-relaxation
      e[0] /= par.relax;
    }
  }
  // Assembles linear system
  // fcu: field from previous iteration [a]
  // ffv: volume flux
  // Output:
  // fcl: linear system
  void Assemble(const FieldCell<Scal>& fcu, const FieldFaceb<Scal>& ffv) {
    Assemble(fcu, ffv, fcucs_);
  }
  void MakeIteration() {
    auto sem = m.GetSem("convdiff-iter");

    auto& prev = fcu_.iter_prev;
    auto& curr = fcu_.iter_curr;

    if (sem("init")) {
      prev.swap(curr);
    }
    if (sem("assemble")) {
      Assemble(prev, *owner_->ffv_, fcucs_);
    }
    if (sem.Nested("solve")) {
      Solve(fcucs_, nullptr, curr, M::LS::T::gen, m);
    }
    if (sem("apply")) {
      // apply, store result in curr
      for (auto c : eb.Cells()) {
        curr[c] += prev[c];
      }
      CalcError();
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
    FieldCell<Scal> fc(eb, 0);
    for (auto c : eb.Cells()) {
      fc[c] = fcucs_[c][0];
    }
    return fc;
  }
  FieldCell<Scal> GetConst() const {
    FieldCell<Scal> fc(eb, 0);
    for (auto c : eb.Cells()) {
      fc[c] = fcucs_[c][Expr::dim - 1];
    }
    return fc;
  }

  Owner* owner_;
  Par par;
  M& m;
  const EB& eb;

  StepData<FieldCell<Scal>> fcu_; // field
  const MapEmbed<BCond<Scal>>& mebc_; // boundary conditions

  FieldCell<Expr> fcucs_; // linear system for correction

  Scal dtp_; // dt prev
  Scal er_; // error
};

template <class EB_>
ConvDiffScalImp<EB_>::ConvDiffScalImp(
    M& m, const EB& eb, const FieldCell<Scal>& fcu,
    const MapEmbed<BCond<Scal>>& mebc, const FieldCell<Scal>* fcr,
    const FieldFaceb<Scal>* ffd, const FieldCell<Scal>* fcs,
    const FieldFaceb<Scal>* ffv, double t, double dt, Par par)
    : Base(t, dt, m, eb, par, fcr, ffd, fcs, ffv)
    , imp(new Imp(this, fcu, mebc)) {}

template <class EB_>
ConvDiffScalImp<EB_>::~ConvDiffScalImp() = default;

template <class EB_>
void ConvDiffScalImp<EB_>::Assemble(
    const FieldCell<Scal>& fcu, const FieldFaceb<Scal>& ffv) {
  imp->Assemble(fcu, ffv);
}

template <class EB_>
void ConvDiffScalImp<EB_>::CorrectField(Step l, const FieldCell<Scal>& uc) {
  imp->CorrectField(l, uc);
}

template <class EB_>
auto ConvDiffScalImp<EB_>::GetDiag() const -> FieldCell<Scal> {
  return imp->GetDiag();
}

template <class EB_>
auto ConvDiffScalImp<EB_>::GetConst() const -> FieldCell<Scal> {
  return imp->GetConst();
}

template <class EB_>
void ConvDiffScalImp<EB_>::StartStep() {
  imp->StartStep();
}

template <class EB_>
void ConvDiffScalImp<EB_>::MakeIteration() {
  imp->MakeIteration();
}

template <class EB_>
void ConvDiffScalImp<EB_>::FinishStep() {
  imp->FinishStep();
}

template <class EB_>
double ConvDiffScalImp<EB_>::GetError() const {
  return imp->GetError();
}

template <class EB_>
auto ConvDiffScalImp<EB_>::GetField(Step l) const -> const FieldCell<Scal>& {
  return imp->fcu_.Get(l);
}
