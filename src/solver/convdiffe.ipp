#include <cmath>
#include <sstream>
#include <stdexcept>

#include "convdiffe.h"
#include "debug/isnan.h"

namespace solver {

template <class M_>
struct ConvDiffScalExp<M_>::Imp {
  using Owner = ConvDiffScalExp<M_>;

  Imp(Owner* owner, const FieldCell<Scal>& fcu,
      const MapFace<std::shared_ptr<CondFace>>& mfc, 
      const MapCell<std::shared_ptr<CondCell>>& mcc,
      std::shared_ptr<Par> par)
      : owner_(owner), par(par), m(owner_->m)
      , mfc_(mfc), mcc_(mcc), dtp_(-1.), er_(0)
  {
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
  void Assemble(const FieldCell<Scal>& fcu, const FieldFace<Scal>& ffv,
                FieldCell<Expr>& fcl) {
    auto sem = m.GetSem("assemble");

    if (sem("assemble")) {
      FieldCell<Vect> fcg = Gradient(Interpolate(fcu, mfc_, m), m);

      FieldFace<Scal> ffq;  // flux tmp
      FieldCell<Scal> fca;  // accumulated

      fca.Reinit(m, 0.);

      // Calc convective fluxes
      Interpolate(fcu, fcg, mfc_, ffv, m, par->sc, par->th, ffq);
      for (auto f : m.Faces()) {
        ffq[f] *= (*owner_->ffv_)[f];
      }
      // Append
      for (IdxCell c : m.Cells()) {
        Scal s = 0.; // sum
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetNeighbourFace(c, q);
          s += ffq[f] * m.GetOutwardFactor(c, q);
        }
        fca[c] += s / m.GetVolume(c) * (*owner_->fcr_)[c];
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
            IdxFace f = m.GetNeighbourFace(c, q);
            s += ffq[f] * m.GetOutwardFactor(c, q);
          }
          fca[c] += s / m.GetVolume(c);
        }
      }

      // time derivative coeffs
      Scal dt = owner_->GetTimeStep();
      std::vector<Scal> ac = GetGradCoeffs(
          0., {-(dt + dtp_), -dt, 0.}, par->second ? 0 : 1);

      // Init system
      fcl.Reinit(m);
      for (IdxCell c : m.Cells()) {
        fca[c] += -(*owner_->fcs_)[c]; // source

        Scal r = (*owner_->fcr_)[c];

        Expr& e = fcl[c];
        e.Clear();
        e.InsertTerm(ac[2] * r, c);
        e.SetConstant(fca[c] +
            (ac[0] * fcu_.time_prev[c] + ac[1] * fcu_.time_curr[c]) * r);
        // Convert to delta-form
        e.SetConstant(e.Evaluate(fcu));
        // Apply under-relaxation
        e[e.Find(c)].a /= par->relax;
      }

      // Overwrite with cell conditions 
      for (auto it = mcc_.cbegin(); it != mcc_.cend(); ++it) {
        IdxCell c(it->GetIdx());
        CondCell* cb = it->GetValue().get(); // cond base
        auto& e = fcl[c];
        if (auto cd = dynamic_cast<CondCellVal<Scal>*>(cb)) {
          Scal v = cd->GetValue() - fcu[c];
          e.SetKnownValueDiag(c, v);
        }
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
    if (sem("copy")) {
      fcu.Reinit(m);
      for (auto c : m.Cells()) {
        fcu[c] = -fcl[c].GetConstant() / fcl[c][0].a;
      }
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
  void CorrectField(Layers l, const FieldCell<Scal>& uc) {
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
  std::shared_ptr<Par> par;
  M& m; // mesh

  LayersData<FieldCell<Scal>> fcu_; // field
  MapFace<std::shared_ptr<CondFace>> mfc_; // face cond
  MapCell<std::shared_ptr<CondCell>> mcc_; // cell cond

  FieldCell<Expr> fcucs_;  // linear system for correction

  Scal dtp_; // dt prev
  Scal er_; // error 
};

template <class M_>
ConvDiffScalExp<M_>::ConvDiffScalExp(
    M& m, const FieldCell<Scal>& fcu, 
    const MapFace<std::shared_ptr<CondFace>>& mfc, 
    const MapCell<std::shared_ptr<CondCell>>& mcc, 
    const FieldCell<Scal>* fcr, const FieldFace<Scal>* ffd,
    const FieldCell<Scal>* fcs, const FieldFace<Scal>* ffv,
    double t, double dt, std::shared_ptr<Par> par)
    : ConvDiffScal<M>(t, dt, m, fcr, ffd, fcs, ffv)
    , imp(new Imp(this, fcu, mfc, mcc, par))
{}

template <class M_>
ConvDiffScalExp<M_>::~ConvDiffScalExp() = default;

template <class M_>
auto ConvDiffScalExp<M_>::GetPar() -> Par* {
  return imp->par.get();
}

template <class M_>
void ConvDiffScalExp<M_>::Assemble(const FieldCell<Scal>& fcu, 
                                   const FieldFace<Scal>& ffv) {
  imp->Assemble(fcu, ffv);
}

template <class M_>
void ConvDiffScalExp<M_>::CorrectField(Layers l, const FieldCell<Scal>& uc) {
  imp->CorrectField(l, uc);
}

template <class M_>
auto ConvDiffScalExp<M_>::GetDiag() const -> FieldCell<Scal> {
  return imp->GetDiag();
}

template <class M_>
auto ConvDiffScalExp<M_>::GetConst() const -> FieldCell<Scal> {
  return imp->GetConst();
}


template <class M_>
void ConvDiffScalExp<M_>::StartStep() {
  imp->StartStep();
}

template <class M_>
void ConvDiffScalExp<M_>::MakeIteration() {
  imp->MakeIteration();
}

template <class M_>
void ConvDiffScalExp<M_>::FinishStep() {
  imp->FinishStep();
}

template <class M_>
double ConvDiffScalExp<M_>::GetError() const {
  return imp->GetError();
}

template <class M_>
auto ConvDiffScalExp<M_>::GetField(Layers l) const -> const FieldCell<Scal>& {
  return imp->fcu_.Get(l);
}


} // namespace solver
