#include <cmath>
#include <sstream>
#include <stdexcept>

#include "convdiffi.h"
#include "debug/isnan.h"

namespace solver {

template <class M_>
struct ConvDiffScalImp<M_>::Imp {
  using Owner = ConvDiffScalImp<M_>;

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
    CHECKNAN(fcu_.time_curr)

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
  // fcl: linear system
  void Assemble(const FieldCell<Scal>& fcu, const FieldFace<Scal>& ffv,
                FieldCell<Expr>& fcl) {
    auto sem = m.GetSem("assemble");

    if (sem("assemble")) {
      FieldCell<Vect> fcg = Gradient(Interpolate(fcu, mfc_, m), m);

      FieldFace<Expr> ffq;  // flux tmp

      // Calc convective fluxes:
			// all inner
      InterpolateI(fcu, fcg, ffv, ffq, m, par->sc, par->df, par->th);
      for (auto f : m.Faces()) {
        ffq[f] *= (*owner_->ffv_)[f];
      }

			// overwrite with bc
      FaceValB<M, Expr> ub(m, mfc_); 
      for (auto it = mfc_.cbegin(); it != mfc_.cend(); ++it) {
        IdxFace f = it->GetIdx();
        Expr e = ub.GetExpr(f);
        ffq[f] = e * (*owner_->ffv_)[f];
      }

      // Init system with convective flux, time derivative, source
      fcl.Reinit(m);
      Scal dt = owner_->GetTimeStep();
      std::vector<Scal> ac = GetGradCoeffs(
          0., {-(dt + dtp_), -dt, 0.}, par->second ? 0 : 1);
      for (IdxCell c : m.Cells()) {
        Expr& e = fcucs_[c];

        Expr sc; // sum convective
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetNeighbourFace(c, q);
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
        FaceGradB<M, Expr> gb(m, mfc_);
        for (auto it = mfc_.cbegin(); it != mfc_.cend(); ++it) {
          IdxFace f = it->GetIdx();
          Expr e = gb.GetExpr(f);
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
          IdxFace f = m.GetNeighbourFace(c, q);
          sd += ffq[f] * m.GetOutwardFactor(c, q);
        }

        auto vol = m.GetVolume(c);
        e += sd / vol;

        // Convert to delta-form
        e.SetConstant(e.Evaluate(fcu));

        // Apply under-relaxation
        e[e.Find(c)].a /= par->relax;
      }

      // Overwrite with cell conditions 
      for (auto it = mcc_.cbegin(); it != mcc_.cend(); ++it) {
        auto dt = owner_->GetTimeStep();
        IdxCell c(it->GetIdx());
        CondCell* cb = it->GetValue().get(); // cond base
        auto& eqn = fcl[c];
        if (auto cd = dynamic_cast<CondCellVal<Scal>*>(cb)) {
          eqn.Clear();
          // TODO: Revise dt coefficient for fixed-value cell condition
          eqn.InsertTerm(1. / dt, c);
          eqn.SetConstant((fcu[c] - cd->GetValue()) / dt);
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
      CHECKNAN(fcucs_)
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
    CHECKNAN(fcu_.time_curr)
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
  const FieldCell<Expr>& GetEquations() const {
    return fcucs_;
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
ConvDiffScalImp<M_>::ConvDiffScalImp(
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
ConvDiffScalImp<M_>::~ConvDiffScalImp() = default;

template <class M_>
auto ConvDiffScalImp<M_>::GetPar() -> Par* {
  return imp->par.get();
}

template <class M_>
void ConvDiffScalImp<M_>::Assemble(const FieldCell<Scal>& fcu, 
                                   const FieldFace<Scal>& ffv) {
  imp->Assemble(fcu, ffv);
}

template <class M_>
void ConvDiffScalImp<M_>::CorrectField(Layers l, const FieldCell<Scal>& uc) {
  imp->CorrectField(l, uc);
}

template <class M_>
auto ConvDiffScalImp<M_>::GetEquations() const -> const FieldCell<Expr>& {
  return imp->GetEquations();
}

template <class M_>
void ConvDiffScalImp<M_>::StartStep() {
  imp->StartStep();
}

template <class M_>
void ConvDiffScalImp<M_>::MakeIteration() {
  imp->MakeIteration();
}

template <class M_>
void ConvDiffScalImp<M_>::FinishStep() {
  imp->FinishStep();
}

template <class M_>
double ConvDiffScalImp<M_>::GetError() const {
  return imp->GetError();
}

template <class M_>
auto ConvDiffScalImp<M_>::GetField(Layers l) const -> const FieldCell<Scal>& {
  return imp->fcu_.Get(l);
}


} // namespace solver
