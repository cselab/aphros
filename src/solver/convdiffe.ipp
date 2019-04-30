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

    // initial guess by extrapolation
    Scal ge = par->guessextra;
    if (ge != 0.) {
      for (auto c : m.AllCells()) {
        fcu_.iter_curr[c] = 
            fcu_.time_curr[c] * (1. + ge) - fcu_.time_prev[c] * ge;
      }
    } else {
      fcu_.iter_curr = fcu_.time_curr;
    }

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
  void Assemble(const FieldCell<Scal>& fcu, const FieldFace<Scal>& ffv,
                FieldCell<Scal>& fcla, FieldCell<Scal>& fclb) {
    auto sem = m.GetSem("assemble");

    if (sem("assemble")) {
      FieldCell<Vect> fcg = Gradient(Interpolate(fcu, mfc_, m), m);

      fcla.Reinit(m);
      fclb.Reinit(m, 0.);

      FieldFace<Scal> ffq;  // flux tmp

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
            IdxFace f = m.GetNeighbourFace(c, q);
            s += ffq[f] * m.GetOutwardFactor(c, q);
          }
          fclb[c] += s / m.GetVolume(c);
        }
      }

      // time derivative coeffs
      Scal dt = owner_->GetTimeStep();
      std::vector<Scal> ac = GetGradCoeffs(
          0., {-(dt + dtp_), -dt, 0.}, par->second ? 0 : 1);

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
        fcla[c] /= par->relax;
      }

      // Overwrite with cell conditions 
      for (auto it = mcc_.cbegin(); it != mcc_.cend(); ++it) {
        IdxCell c(it->GetIdx());
        CondCell* cb = it->GetValue().get(); // cond base
        if (auto cd = dynamic_cast<CondCellVal<Scal>*>(cb)) {
          fcla[c] = 1.;
          fclb[c] = cd->GetValue() - fcu[c];
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
    Assemble(fcu, ffv, fcla_, fclb_);
  }
  // Solves linear system.
  // fcla * u + fclb = 0: system to solve
  // Output:
  // fcu: result
  void Solve(const FieldCell<Scal>& fcla, const FieldCell<Scal>& fclb,
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
    return fcla_;
  }
  FieldCell<Scal> GetConst() const {
    return fclb_;
  }

  Owner* owner_;
  std::shared_ptr<Par> par;
  M& m; // mesh

  LayersData<FieldCell<Scal>> fcu_; // field
  MapFace<std::shared_ptr<CondFace>> mfc_; // face cond
  MapCell<std::shared_ptr<CondCell>> mcc_; // cell cond

  // diagonal linear system: fcla * (fcup - fcu) + fclb = 0
  FieldCell<Scal> fcla_;
  FieldCell<Scal> fclb_;

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
