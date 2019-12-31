#include <cmath>
#include <sstream>
#include <stdexcept>

#include "approx.h"
#include "convdiffe_eb.h"
#include "debug/isnan.h"

template <class M_>
struct ConvDiffScalExpEmbed<M_>::Imp {
  using Owner = ConvDiffScalExpEmbed<M_>;
  using Vect = typename M::Vect;
  using EB = Embed<M>;
  using Type = typename EB::Type;

  Imp(Owner* owner, const FieldCell<Scal>& fcu, const MapCondFace& mfc,
      const Embed<M>& eb, const FieldEmbed<Scal>* fed, size_t bc, Scal bcu)
      : owner_(owner)
      , par(owner_->GetPar())
      , m(owner_->m)
      , mfc_(mfc)
      , eb(eb)
      , fed_(fed)
      , bc_(bc)
      , bcu_(bcu)
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
      const FieldCell<Scal>& fcu, const FieldEmbed<Scal>& fev,
      FieldCell<Scal>& fcla, FieldCell<Scal>& fclb) {
    auto sem = m.GetSem("assemble");

    if (sem("assemble")) {
      // FieldCell<Vect> fcg = eb.Gradient(eb.Interpolate(fcu, bc_, bcu_), m);

      fcla.Reinit(m);
      fclb.Reinit(m, 0.);

      // Convective fluxes
      {
        FieldEmbed<Scal> feq = eb.Interpolate(fcu, bc_, bcu_);
        for (auto f : m.Faces()) {
          if (eb.GetType(f) != Type::excluded) {
            feq[f] *= fev[f];
          }
        }
        for (auto c : m.Cells()) {
          if (eb.GetType(c) == Type::cut) {
            feq[c] *= fev[c];
          }
        }
        // Append
        for (IdxCell c : m.Cells()) {
          Scal s = 0.; // sum
          for (auto q : m.Nci(c)) {
            IdxFace f = m.GetFace(c, q);
            s += feq[f] * m.GetOutwardFactor(c, q);
          }
          fclb[c] += s / m.GetVolume(c) * (*owner_->fcr_)[c];
        }
      }

      // Diffusive fluxes
      if (fed_) {
        FieldEmbed<Scal> feq = eb.Gradient(fcu, bc_, bcu_);
        for (auto f : m.Faces()) {
          if (eb.GetType(f) != Type::excluded) {
            feq[f] *= (-(*fed_)[f]) * eb.GetArea(f);
          }
        }
        for (auto c : m.Cells()) {
          if (eb.GetType(c) == Type::cut) {
            feq[c] *= (-(*fed_)[c]) * eb.GetArea(c);
          }
        }
        // Append
        for (IdxCell c : m.Cells()) {
          Scal s = 0.; // sum
          for (auto q : m.Nci(c)) {
            IdxFace f = m.GetFace(c, q);
            s += feq[f] * m.GetOutwardFactor(c, q);
          }
          fclb[c] += s / m.GetVolume(c);
        }
      }

      // time derivative coeffs
      const Scal dt = owner_->GetTimeStep();
      const std::vector<Scal> ac =
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

  void Assemble(
      const FieldCell<Scal>& fcu, const FieldFace<Scal>& ffv,
      FieldCell<Scal>& fcla, FieldCell<Scal>& fclb) {
    const FieldCell<Scal> fcv(m, 0);
    const FieldEmbed<Scal> fev(fcv, ffv);
    Assemble(fcu, fev, fcla, fclb);
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

  StepData<FieldCell<Scal>> fcu_; // field
  const MapCondFace& mfc_; // face cond

  const Embed<M>& eb;
  const FieldEmbed<Scal>* fed_;

  size_t bc_; // boundary conditions, 0: value, 1: gradient
  Scal bcu_; // value or grad.dot.outer_normal

  // diagonal linear system: fcla * (fcup - fcu) + fclb = 0
  FieldCell<Scal> fcla_;
  FieldCell<Scal> fclb_;

  Scal dtp_; // dt prev
  Scal er_; // error
};

template <class M_>
ConvDiffScalExpEmbed<M_>::ConvDiffScalExpEmbed(
    M& m, const Embed<M>& eb, const FieldCell<Scal>& fcu,
    const MapCondFace& mfc, size_t bc, Scal bcu, const FieldCell<Scal>* fcr,
    const FieldEmbed<Scal>* fed, const FieldCell<Scal>* fcs,
    const FieldFace<Scal>* ffv, double t, double dt, Par par)
    : ConvDiffScal<M>(t, dt, m, par, fcr, nullptr, fcs, ffv)
    , imp(new Imp(this, fcu, mfc, eb, fed, bc, bcu)) {}

template <class M_>
ConvDiffScalExpEmbed<M_>::~ConvDiffScalExpEmbed() = default;

template <class M_>
void ConvDiffScalExpEmbed<M_>::Assemble(
    const FieldCell<Scal>& fcu, const FieldFace<Scal>& ffv) {
  imp->Assemble(fcu, ffv);
}

template <class M_>
void ConvDiffScalExpEmbed<M_>::CorrectField(Step l, const FieldCell<Scal>& uc) {
  imp->CorrectField(l, uc);
}

template <class M_>
auto ConvDiffScalExpEmbed<M_>::GetDiag() const -> FieldCell<Scal> {
  return imp->GetDiag();
}

template <class M_>
auto ConvDiffScalExpEmbed<M_>::GetConst() const -> FieldCell<Scal> {
  return imp->GetConst();
}

template <class M_>
void ConvDiffScalExpEmbed<M_>::StartStep() {
  imp->StartStep();
}

template <class M_>
void ConvDiffScalExpEmbed<M_>::MakeIteration() {
  imp->MakeIteration();
}

template <class M_>
void ConvDiffScalExpEmbed<M_>::FinishStep() {
  imp->FinishStep();
}

template <class M_>
double ConvDiffScalExpEmbed<M_>::GetError() const {
  return imp->GetError();
}

template <class M_>
auto ConvDiffScalExpEmbed<M_>::GetField(Step l) const
    -> const FieldCell<Scal>& {
  return imp->fcu_.Get(l);
}
