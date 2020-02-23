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
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;
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
  // Assembles linear system
  // fcu: field from previous iteration [a]
  // ffv: volume flux
  // Output:
  // fcl: linear system
  void Assemble(
      const FieldCell<Scal>& fcu, const FieldFace<Scal>& ffv,
      FieldCell<Expr>& fcl) {
    auto sem = m.GetSem("assemble");

    if (sem("assemble")) {
      const FieldCell<Vect> fcg = Gradient(UEB::Interpolate(fcu, mebc_, m), m);

      FieldFace<Expr> ffq; // flux tmp

      // Calc convective fluxes:
      // all inner
      InterpolateI(fcu, fcg, ffv, ffq, m, par.sc, par.df, par.th);
      for (auto f : m.Faces()) {
        ffq[f] *= (*owner_->ffv_)[f];
      }

      // overwrite with bc
      for (auto& p : mebc_.GetMapFace()) {
        const IdxFace f = p.first;
        Expr e = UEB::template GetExprVal<Expr>(f, p.second, m);
        e.SortTerms();
        ffq[f] = e * (*owner_->ffv_)[f];
      }

      // Init system with convective flux, time derivative, source
      fcl.Reinit(m);
      Scal dt = owner_->GetTimeStep();
      std::vector<Scal> ac =
          GetGradCoeffs(0., {-(dt + dtp_), -dt, 0.}, par.second ? 0 : 1);
      for (IdxCell c : m.Cells()) {
        Expr& e = fcl[c];

        Expr sc; // sum convective
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetFace(c, q);
          sc += ffq[f] * m.GetOutwardFactor(c, q);
        }

        Expr tt; // time derivative term
        tt.InsertTerm(ac[2], c);
        tt.SetConstant(ac[0] * fcu_.time_prev[c] + ac[1] * fcu_.time_curr[c]);

        auto vol = m.GetVolume(c);
        e = (tt + sc / vol) * (*owner_->fcr_)[c] - Expr((*owner_->fcs_)[c]);
      }

      if (owner_->ffd_) {
        // Calc diffusive fluxes
        // all inner
        GradientI(ffq, m);
        for (auto f : m.Faces()) {
          ffq[f] *= (-(*owner_->ffd_)[f]) * m.GetArea(f);
        }
        // overwrite with bc
        for (auto& p : mebc_.GetMapFace()) {
          const IdxFace f = p.first;
          Expr e = UEB::template GetExprGrad<Expr>(f, p.second, m);
          e.SortTerms();
          ffq[f] = e * (-(*owner_->ffd_)[f]) * m.GetArea(f);
        }
      } else {
        ffq.Reinit(m);
      }

      // Append diffusive flux, convert to delta-form, apply underelaxation
      for (IdxCell c : m.Cells()) {
        Expr& e = fcl[c];

        Expr sd; // sum diffusive
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetFace(c, q);
          sd += ffq[f] * m.GetOutwardFactor(c, q);
        }

        auto vol = m.GetVolume(c);
        e += sd / vol;

        // Convert to delta-form
        e.SetConstant(e.Evaluate(fcu));

        // Apply under-relaxation
        e[e.Find(c)].a /= par.relax;
      }
    }
  }
  // Assembles linear system
  // fcu: field from previous iteration [a]
  // ffv: volume flux
  // Output:
  // fcl: linear system
  void Assemble(const FieldCell<Scal>& fcu, const FieldFace<Scal>& ffv) {
    Assemble(fcu, ffv, fcucs_);
  }
  // Solves linear system.
  // fcl : system to solve
  // Output:
  // fcu: result
  // m.GetSolveTmp(): modified temporary fields
  void Solve(const FieldCell<Expr>& fcl, FieldCell<Scal>& fcu) {
    auto sem = m.GetSem("solve");
    if (sem("convert")) {
      std::vector<Scal>* lsa;
      std::vector<Scal>* lsb;
      std::vector<Scal>* lsx;
      m.GetSolveTmp(lsa, lsb, lsx);
      auto l = ConvertLs(fcl, *lsa, *lsb, *lsx, m);
      using T = typename M::LS::T;
      l.t = T::gen;
      m.Solve(l);
    }
    if (sem("copy")) {
      std::vector<Scal>* lsa;
      std::vector<Scal>* lsb;
      std::vector<Scal>* lsx;
      m.GetSolveTmp(lsa, lsb, lsx);

      fcu.Reinit(m);
      size_t i = 0;
      for (auto c : m.Cells()) {
        fcu[c] = (*lsx)[i++];
      }
      assert(i == lsx->size());
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
      Assemble(prev, *owner_->ffv_, fcucs_);
    }
    if (sem.Nested("solve")) {
      Solve(fcucs_, curr); // solve for correction, store in curr
    }
    if (sem("apply")) {
      CHECKNAN(fcucs_, m.CN())
      // calc error (norm of correction)
      Scal er = 0;
      for (auto c : m.Cells()) {
        er = std::max<Scal>(er, curr[c]);
      }
      er_ = er;

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
      m.Comm(&u);
    }
  }
  const FieldCell<Expr>& GetEquations() const {
    return fcucs_;
  }
  FieldCell<Scal> GetDiag() const {
    FieldCell<Scal> fc(m);
    for (auto c : m.Cells()) {
      fc[c] = fcucs_[c].Coeff(c);
    }
    return fc;
  }
  FieldCell<Scal> GetConst() const {
    FieldCell<Scal> fc(m);
    for (auto c : m.Cells()) {
      fc[c] = fcucs_[c].GetConstant();
    }
    return fc;
  }

  Owner* owner_;
  Par par;
  M& m; // mesh
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
    const FieldFace<Scal>* ffd, const FieldCell<Scal>* fcs,
    const FieldFace<Scal>* ffv, double t, double dt, Par par)
    : ConvDiffScal<M>(t, dt, m, eb, par, fcr, ffd, fcs, ffv)
    , imp(new Imp(this, fcu, mebc)) {}

template <class EB_>
ConvDiffScalImp<EB_>::~ConvDiffScalImp() = default;

template <class EB_>
void ConvDiffScalImp<EB_>::Assemble(
    const FieldCell<Scal>& fcu, const FieldFace<Scal>& ffv) {
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
