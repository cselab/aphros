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
  using ExprLegacy = Expression<Scal, IdxCell, 1 + dim * 2>;
  // Expression on face: v[0] * cm + v[1] * cp + v[2]
  using ExprFace = generic::Vect<Scal, 3>;
  // Expression on cell: v[0] * c + v[1] * cxm + ... + v[6] * czp + v[7]
  using Expr = generic::Vect<Scal, M::dim * 2 + 2>;
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
  Scal Eval(const Expr& e, IdxCell c, const FieldCell<Scal>& fcu) {
    Scal r = e[Expr::dim - 1];
    r += fcu[c] * e[0];
    for (auto q : m.Nci(c)) {
      r += fcu[m.GetCell(c, q)] * e[1 + q];
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
    FieldFace<ExprFace> ffq(m, ExprFace(0));
    // convective fluxes
    {
      const FieldCell<Vect> fcg = Gradient(UEB::Interpolate(fcu, mebc_, m), m);
      const FieldFace<ExprFace> ffu = UEB::InterpolateUpwindImplicit(
          fcu, mebc_, par.sc, par.df, fcg, ffv, m);
      for (auto f : eb.Faces()) {
        ffq[f] += ffu[f] * ffv[f];
      }
    }
    // diffusive fluxes
    if (owner_->ffd_) {
      const FieldFace<ExprFace> ffg = UEB::GradientImplicit(fcu, mebc_, m);
      for (auto f : m.Faces()) {
        ffq[f] -= ffg[f] * (*owner_->ffd_)[f] * m.GetArea(f);
      }
    }

    const Scal dt = owner_->GetTimeStep();
    std::vector<Scal> time_coeff =
        GetGradCoeffs(0., {-(dt + dtp_), -dt, 0.}, par.second ? 0 : 1);

    fcl.Reinit(m, Expr(0));
    for (auto c : m.Cells()) {
      Expr sum(0);
      for (auto q : m.Nci(c)) {
        const IdxFace f = m.GetFace(c, q);
        const ExprFace v = ffq[f] * m.GetOutwardFactor(c, q);
        sum[0] += v[1 - q % 2];
        sum[1 + q] += v[q % 2];
        sum[Expr::dim - 1] += v[2];
      }

      Expr td(0); // time derivative
      td[0] = time_coeff[2];
      td[Expr::dim - 1] =
          time_coeff[0] * fcu_.time_prev[c] + time_coeff[1] * fcu_.time_curr[c];

      Expr& e = fcl[c];
      e = (td + sum / m.GetVolume(c)) * (*owner_->fcr_)[c];
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
  // Solves linear system.
  // fcl : system to solve
  // Output:
  // fcu: result
  // m.GetSolveTmp(): modified temporary fields
  void Solve(const FieldCell<ExprLegacy>& fcl, FieldCell<Scal>& fcu) {
    auto sem = m.GetSem("solve");
    if (sem("convert")) {
      std::vector<Scal>* lsa;
      std::vector<Scal>* lsb;
      std::vector<Scal>* lsx;
      m.GetSolveTmp(lsa, lsb, lsx);
      auto l = ConvertLs(fcl, *lsa, *lsb, *lsx, m);
      l.t = M::LS::T::gen;
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
  // Solve linear system fce = 0
  // fce: expressions [i]
  // Output:
  // fc: result [a]
  // m.GetSolveTmp(): modified temporary fields
  void Solve(const FieldCell<Expr>& fce, FieldCell<Scal>& fc) {
    auto sem = m.GetSem("solve");
    if (sem("solve")) {
      std::vector<Scal>* lsa;
      std::vector<Scal>* lsb;
      std::vector<Scal>* lsx;
      m.GetSolveTmp(lsa, lsb, lsx);
      lsx->resize(m.GetInBlockCells().size());
      size_t i = 0;
      for (auto c : m.Cells()) {
        (void)c;
        (*lsx)[i++] = 0;
      }
      auto l = ConvertLsCompact(fce, *lsa, *lsb, *lsx, m);
      l.t = M::LS::T::gen;
      m.Solve(l);
    }
    if (sem("copy")) {
      std::vector<Scal>* lsa;
      std::vector<Scal>* lsb;
      std::vector<Scal>* lsx;
      m.GetSolveTmp(lsa, lsb, lsx);

      fc.Reinit(m);
      size_t i = 0;
      for (auto c : m.Cells()) {
        fc[c] = (*lsx)[i++];
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
    if (sem("assemble")) {
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
  FieldCell<Scal> GetDiag() const {
    FieldCell<Scal> fc(m);
    for (auto c : m.Cells()) {
      fc[c] = fcucs_[c][0];
    }
    return fc;
  }
  FieldCell<Scal> GetConst() const {
    FieldCell<Scal> fc(m);
    for (auto c : m.Cells()) {
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
