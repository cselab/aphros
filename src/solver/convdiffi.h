#pragma once

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <memory>

#include "convdiff.h"

namespace solver {

template <class M_>
class ConvectionDiffusionScalarImplicit : public ConvectionDiffusionScalar<M_> {
 public:
  using M = M_; // mesh
  using P = ConvectionDiffusionScalar<M>; // parent
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;

  struct Par {
    Scal relax = 1.;      // relaxation factor [0,1] (1 -- no relaxation)
    Scal guessextra = 0.; // next iteration guess extrapolation weight [0,1]
    bool second = true; // second order in time
    ConvSc sc = ConvSc::quick; // scheme for convective flux (see convdiffi.h)
    Scal df = 1.; // deferred correction factor
    Scal th = 1e-10; // threshold for flow direction
  };

  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;
  ConvectionDiffusionScalarImplicit(
      M& m,
      const FieldCell<Scal>& fcu, // initial field
      const MapFace<std::shared_ptr<CondFace>>& mfc, // face conditions
      const MapCell<std::shared_ptr<CondCell>>& mcc, // cell conditions
      const FieldCell<Scal>* fcr, // density
      const FieldFace<Scal>* ffd, // diffusion
      const FieldCell<Scal>* fcs, // source
      const FieldFace<Scal>* ffv, // volume flux
      double t, double dt, std::shared_ptr<Par> par)
      : ConvectionDiffusionScalar<M>(t, dt, fcr, ffd, fcs, ffv)
      , m(m), par(par), mfc_(mfc), mcc_(mcc), dtp_(-1.), er_(0)
  {
    fcu_.time_curr = fcu;
    fcu_.time_prev = fcu;
    fcu_.iter_curr = fcu;
  }
  Par* GetPar() { return par.get(); }
  // Fields:
  void StartStep() override {
    this->ClearIter();
    CheckNan(fcu_.time_curr, "fcu_.time_curr", m);

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
      dtp_ = this->GetTimeStep();
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
        ffq[f] *= (*ffv_)[f];
      }

			// overwrite with bc
      FaceValB<M, Expr> ub(m, mfc_); 
      for (auto it = mfc_.cbegin(); it != mfc_.cend(); ++it) {
        IdxFace f = it->GetIdx();
        Expr e = ub.GetExpr(f);
        ffq[f] = e * (*ffv_)[f];
      }

      // Init system with convective flux, time derivative, source
      fcl.Reinit(m);
      Scal dt = this->GetTimeStep();
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
        e = (tt + sc / vol) * (*fcr_)[c] - Expr((*fcs_)[c]);
      }

      // Calc diffusive fluxes
      // all inner
      GradientI(ffq, m);
      for (auto f : m.Faces()) {
        ffq[f] *= (-(*ffd_)[f]) * m.GetArea(f);
      }
			// overwrite with bc
      FaceGradB<M, Expr> gb(m, mfc_);
      for (auto it = mfc_.cbegin(); it != mfc_.cend(); ++it) {
        IdxFace f = it->GetIdx();
				Expr e = gb.GetExpr(f);
        ffq[f] = e * (-(*ffd_)[f]) * m.GetArea(f);
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
        e[e.Find(c)].coeff /= par->relax;
      }

      // Overwrite with cell conditions 
      for (auto it = mcc_.cbegin(); it != mcc_.cend(); ++it) {
        auto dt = this->GetTimeStep();
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
  // lsa_, lsb_, lsx_: modified temporary fields
  void Solve(const FieldCell<Expr>& fcl, FieldCell<Scal>& fcu) {
    auto sem = m.GetSem("solve");
    if (sem("convert")) {
      auto l = ConvertLs(fcl, lsa_, lsb_, lsx_, m);
      using T = typename M::LS::T;
      l.t = T::gen;
      m.Solve(l);
    }
    if (sem("copy")) {
      fcu.Reinit(m);
      size_t i = 0;
      for (auto c : m.Cells()) {
        fcu[c] = lsx_[i++];
      }
      assert(i == lsx_.size());
    }
  }
  void MakeIteration() override {
    auto sem = m.GetSem("convdiff-iter");

    auto& prev = fcu_.iter_prev;
    auto& curr = fcu_.iter_curr;

    if (sem("init")) {
      prev.swap(curr);
    }
    if (sem.Nested("assemble")) {
      Assemble(prev, *ffv_, fcucs_);
    }
    if (sem.Nested("solve")) {
      Solve(fcucs_, curr); // solve for correction, store in curr
    }
    if (sem("apply")) {
      CheckNan(curr, "convdiffi:fcuc_", m);
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
      this->IncIter();
    }
  }
  void FinishStep() override {
    fcu_.time_prev.swap(fcu_.time_curr);
    fcu_.time_curr = fcu_.iter_curr;
    CheckNan(fcu_.time_curr, "fcu_.time_curr", m);
    this->IncTime();
    dtp_ = this->GetTimeStep();
  }
  double GetError() const override {
    return er_;
  }
  const FieldCell<Scal>& GetField() const override {
    return fcu_.time_curr;
  }
  const FieldCell<Scal>& GetField(Layers layer) const override {
    return fcu_.Get(layer);
  }
  // Apply correction to field and comm
  // uc: correction [i]
  // Output:
  // u(l) += uc [a]
  void CorrectField(Layers l, const FieldCell<Scal>& uc) override {
    auto sem = m.GetSem("corr");
    if (sem("apply")) {
      auto& u = fcu_.Get(l);
      for (auto c : m.Cells()) {
        u[c] += uc[c];
      }
      m.Comm(&u);
    }
  }
  const FieldCell<Expr>& GetEquations() const override {
    return fcucs_;
  }

 private:
  M& m; // mesh
  std::shared_ptr<Par> par; // parameters
  LayersData<FieldCell<Scal>> fcu_; // field
  MapFace<std::shared_ptr<CondFace>> mfc_; // face cond
  MapCell<std::shared_ptr<CondCell>> mcc_; // cell cond

  using P::fcr_;
  using P::ffd_;
  using P::fcs_;
  using P::ffv_;

  FieldCell<Expr> fcucs_;  // linear system for correction

  // tmp for ConvertLs
  std::vector<Scal> lsa_, lsb_, lsx_;

  Scal dtp_; // dt prev
  Scal er_; // error 
};

} // namespace solver
