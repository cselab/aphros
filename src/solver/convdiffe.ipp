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

  Imp(Owner* owner, const Args& args)
      : owner_(owner)
      , par(owner_->GetPar())
      , m(owner_->m)
      , eb(owner_->eb)
      , mebc_(args.mebc)
      , dtp_(-1.)
      , error_(0) {
    fcu_.time_curr = args.fcu;
    fcu_.time_prev = args.fcu;
    fcu_.iter_curr = args.fcu;
  }
  // Fields:
  void StartStep() {
    owner_->ClearIter();
    CHECKNAN(fcu_.time_curr, m.flags.check_nan)

    fcu_.iter_curr = fcu_.time_curr;

    if (dtp_ == -1.) { // TODO: revise
      dtp_ = owner_->GetTimeStep();
    }

    error_ = 0.;
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
    fcla.Reinit(m, 0);
    fclb.Reinit(m, 0);
    // convective fluxes
    if (!par.stokes) {
      FieldFaceb<Scal> ffq(eb, 0);
      const FieldCell<Vect> fcg =
          UEB::Gradient(UEB::Interpolate(fcu, mebc_, eb), eb);
      const FieldFaceb<Scal> ffu =
          UEB::InterpolateUpwind(fcu, mebc_, par.sc, fcg, ffv, eb);
      eb.LoopFaces([&](auto cf) { //
        ffq[cf] = ffu[cf] * ffv[cf];
      });
      for (auto c : eb.Cells()) {
        Scal sum = 0;
        eb.LoopNci(c, [&](auto q) {
          const auto cf = eb.GetFace(c, q);
          sum += ffq[cf] * eb.GetOutwardFactor(c, q);
        });
        fclb[c] += sum * (*owner_->fcr_)[c];
      }
    }
    // diffusive fluxes
    if (owner_->ffd_) {
      FieldFaceb<Scal> ffq(eb, 0);
      const FieldFaceb<Scal> ffg = UEB::Gradient(fcu, mebc_, eb);
      eb.LoopFaces([&](auto cf) { //
        ffq[cf] = -ffg[cf] * (*owner_->ffd_)[cf] * eb.GetArea(cf);
      });
      for (auto c : eb.Cells()) {
        Scal sum = 0;
        eb.LoopNci(c, [&](auto q) {
          const auto cf = eb.GetFace(c, q);
          sum += ffq[cf] * eb.GetOutwardFactor(c, q);
        });
        fclb[c] += sum;
      }
    }
    fclb = UEB::RedistributeCutCells(fclb, eb);

    // time derivative
    if (!par.stokes) {
      const Scal dt = owner_->GetTimeStep();
      const std::vector<Scal> tt =
          GetGradCoeffs(0., {-(dt + dtp_), -dt, 0.}, par.second ? 0 : 1);
      for (auto c : eb.Cells()) {
        const Scal a = eb.GetVolume(c) * (*owner_->fcr_)[c];
        fcla[c] += tt[2] * a;
        fclb[c] += (tt[0] * fcu_.time_prev[c] + tt[1] * fcu_.time_curr[c]) * a;
      }
    } else {
      for (auto c : eb.Cells()) {
        fcla[c] += eb.GetVolume(c);
        fclb[c] -= fcu[c] * eb.GetVolume(c);
      }
    }

    for (auto c : eb.Cells()) {
      // source
      fclb[c] -= (*owner_->fcs_)[c] * eb.GetVolume(c);
      // delta form
      fclb[c] += fcla[c] * fcu[c];
      // under-relaxation
      fcla[c] /= par.relax;
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
    if (sem("assemble")) {
      Assemble(prev, *owner_->ffv_, fcla_, fclb_);
    }
    if (sem("solve")) {
      Solve(fcla_, fclb_, curr); // solve for correction, store in curr
    }
    if (sem("apply")) {
      // apply, store result in curr
      for (auto c : eb.Cells()) {
        curr[c] += prev[c];
      }
      error_ = CalcError();
      m.Reduce(&error_, "max");
      m.Comm(&curr);
      owner_->IncIter();
    }
  }
  void FinishStep() {
    fcu_.time_prev.swap(fcu_.time_curr);
    fcu_.time_curr = fcu_.iter_curr;
    CHECKNAN(fcu_.time_curr, m.flags.check_nan)
    owner_->IncTime();
    dtp_ = owner_->GetTimeStep();
  }
  Scal CalcError() {
    // max difference between iter_curr and iter_prev
    auto& prev = fcu_.iter_prev;
    auto& curr = fcu_.iter_curr;
    Scal a = 0;
    for (auto c : eb.Cells()) {
      a = std::max<Scal>(a, std::abs(curr[c] - prev[c]));
    }
    return a;
  }
  double GetError() const {
    return error_;
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
      fc[c] = fcla_[c] / eb.GetVolume(c);
    }
    return fc;
  }
  FieldCell<Scal> GetConst() const {
    FieldCell<Scal> fc(eb, 0);
    for (auto c : eb.Cells()) {
      fc[c] = fclb_[c] / eb.GetVolume(c);
    }
    return fc;
  }

  Owner* owner_;
  Par par;
  M& m; // mesh
  const EB& eb; // embed mesh

  StepData<FieldCell<Scal>> fcu_; // field
  const MapEmbed<BCond<Scal>>& mebc_; // boundary conditions

  // diagonal linear system: fcla * (fcup - fcu) + fclb = 0
  FieldCell<Scal> fcla_;
  FieldCell<Scal> fclb_;

  Scal dtp_; // dt prev
  Scal error_; // error
};

template <class EB_>
ConvDiffScalExp<EB_>::ConvDiffScalExp(M& m_, const EB& eb_, const Args& args)
    : Base(
          args.t, args.dt, m_, eb_, args.par, args.fcr, args.ffd, args.fcs,
          args.ffv)
    , imp(new Imp(this, args)) {}

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
