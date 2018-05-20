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
      , m(m), par(par), mfc_(mfc), mcc_(mcc), dtp_(-1.)
  {
    fcu_.time_curr = fcu;
    fcu_.time_prev = fcu_.time_curr;
  }
  Par* GetPar() { return par.get(); }
  void StartStep() override {
    this->GetIter();
    if (IsNan(fcu_.time_curr)) {
      throw std::runtime_error("NaN initial field");
    }
    fcu_.iter_curr = fcu_.time_curr;
    Scal ge = par->guessextra;
    for (auto c : m.Cells()) {
      fcu_.iter_curr[c] +=
          (fcu_.time_curr[c] - fcu_.time_prev[c]) * ge;
    }
    if (dtp_ == -1.) {
      dtp_ = this->GetTimeStep();
    }
  }
  // Assembles linear system
  // fcu: field from previous iteration [a]
  // ffv: volume flux
  // Output:
  // fcucs_: linear system
  void Assemble(const FieldCell<Scal>& fcu, const FieldFace<Scal>& ffv) {
    auto sem = m.GetSem("assemble");

    if (sem("assemble")) {
      fcgu_ = Gradient(Interpolate(fcu, mfc_, m), m);

      // Calc convective fluxes:
			// all inner
      InterpolateI(fcu, fcgu_, ffv, ffqc_, m, par->sc);
      for (auto f : m.Faces()) {
        ffqc_[f] *= (*ffv_)[f];
      }

			// overwrite with bc
      FaceValB<M, Expr> ub(m, mfc_); 
      for (auto it = mfc_.cbegin(); it != mfc_.cend(); ++it) {
        IdxFace f = it->GetIdx();
        Expr e = ub.GetExpr(f);
        ffqc_[f] = e * (*ffv_)[f];
      }

      // Calc diffusive fluxes
      // all inner
      GradientI(ffqd_, m);
      for (auto f : m.Faces()) {
        ffqd_[f] *= (-(*ffd_)[f]) * m.GetArea(f);
      }
			// overwrite with bc
      FaceGradB<M, Expr> gb(m, mfc_);
      for (auto it = mfc_.cbegin(); it != mfc_.cend(); ++it) {
        IdxFace f = it->GetIdx();
				Expr e = gb.GetExpr(f);
        ffqd_[f] = e * (-(*ffd_)[f]) * m.GetArea(f);
      }

      // Assemble the system
      fcucs_.Reinit(m);
      for (IdxCell c : m.Cells()) {
        Expr& e = fcucs_[c];
        Expr sc; // sum convective
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetNeighbourFace(c, q);
          sc += ffqc_[f] * m.GetOutwardFactor(c, q);
        }

        Expr sd; // sum diffusive
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetNeighbourFace(c, q);
          sd += ffqd_[f] * m.GetOutwardFactor(c, q);
        }

        Scal dt = this->GetTimeStep();
        std::vector<Scal> ac = GetGradCoeffs(
            0., {-(dt + dtp_), -dt, 0.}, par->second ? 0 : 1);

        Expr tt; // time derivative term
        tt.InsertTerm(ac[2], c);
        tt.SetConstant(ac[0] * fcu_.time_prev[c] + ac[1] * fcu_.time_curr[c]);

        auto vol = m.GetVolume(c);
        e = (tt + sc / vol) * (*fcr_)[c] + sd / vol - Expr((*fcs_)[c]);

        // Convert to delta-form
        e.SetConstant(e.Evaluate(fcu));

        // Apply under-relaxation
        e[e.Find(c)].coeff /= par->relax;
      }

      // Include cell conditions for velocity
      for (auto it = mcc_.cbegin(); it != mcc_.cend(); ++it) {
        auto dt = this->GetTimeStep();
        IdxCell c(it->GetIdx());
        CondCell* cb = it->GetValue().get(); // cond base
        auto& eqn = fcucs_[c];
        if (auto cd = dynamic_cast<CondCellVal<Scal>*>(cb)) {
          eqn.Clear();
          // TODO: Revise dt coefficient for fixed-value cell condition
          eqn.InsertTerm(1. / dt, c);
          eqn.SetConstant((fcu[c] - cd->GetValue()) / dt);
        }
      }
    }

  }
  // Solve linear system.
  // fcucs_ : system to solve
  // Output:
  // fcu: result
  // lsa_, lsb_, lsx_: modified temporary fields
  void Solve(FieldCell<Scal>& fcu) {
    auto sem = m.GetSem("solve");
    if (sem("convert")) {
      auto l = ConvertLs(fcucs_, lsa_, lsb_, lsx_, m);
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
    if (sem.Nested("init")) {
      fcu_.iter_prev = fcu_.iter_curr;
    }
    if (sem.Nested("assemble")) {
      Assemble(fcu_.iter_prev, *ffv_);
    }
    if (sem.Nested("solve")) {
      Solve(fcuc_);
    }
    if (sem("apply")) {
      auto& fc_prev = fcu_.iter_prev;
      auto& fc_curr = fcu_.iter_curr;
      for (auto c : m.Cells()) {
        fc_curr[c] = fc_prev[c] + fcuc_[c];
      }
      m.Comm(&fc_curr);
      this->IncIter();
    }
  }
  void FinishStep() override {
    fcu_.time_prev = fcu_.time_curr;
    fcu_.time_curr = fcu_.iter_curr;
    if (IsNan(fcu_.time_curr)) {
      throw std::runtime_error("NaN field");
    }
    this->IncTime();
    dtp_ = this->GetTimeStep();
  }
  double GetError() const override {
    if (this->GetIter() == 0) {
      return 1.;
    }
    return CalcDiff(fcu_.iter_curr, fcu_.iter_prev, m);
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
  // u(l) += du [a]
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

  FieldFace<Expr> ffqc_;   // convective flux
  FieldFace<Expr> ffqd_;   // diffusive flux
  FieldCell<Expr> fcucs_;  // linear system for correction
  FieldCell<Scal> fcuc_;   // field correction
  FieldCell<Vect> fcgu_;   // gradient of field

  // tmp for ConvertLs
  std::vector<Scal> lsa_, lsb_, lsx_;

  Scal dtp_; // dt prev
};

} // namespace solver
