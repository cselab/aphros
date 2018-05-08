#pragma once

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <memory>

#include "convdiff.h"

namespace solver {

template <class M_>
class ConvectionDiffusionScalarImplicit :
    public ConvectionDiffusionScalar<M_> {
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;

  M& m;
  LayersData<FieldCell<Scal>> fc_field_;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

  MapFace<std::shared_ptr<ConditionFace>> mfc_;
  MapCell<std::shared_ptr<ConditionCell>> mcc_;

  // Common buffers:
  FieldFace<Expr> ff_cflux_;
  FieldFace<Expr> ff_dflux_;
  FieldCell<Expr> fc_system_;
  FieldCell<Scal> fc_corr_;
  FieldCell<Vect> fc_grad_;
  // LS
  std::vector<Scal> lsa_, lsb_, lsx_;

 public:
  struct Par {
    Scal relax = 1.;      // relaxation factor [0,1] (1 -- no relaxation)
    Scal guessextra = 0.; // next iteration guess extrapolation weight [0,1]
    bool second = true; // second order in time
  };
  std::shared_ptr<Par> par;
  Par* GetPar() { return par.get(); }
  ConvectionDiffusionScalarImplicit(
      M& m,
      const FieldCell<Scal>& fcu, // initial field
      const MapFace<std::shared_ptr<ConditionFace>>& mfc, // face conditions
      const MapCell<std::shared_ptr<ConditionCell>>& mcc, // cell conditions
      const FieldCell<Scal>* fcr, // density
      const FieldFace<Scal>* ffd, // diffusion
      const FieldCell<Scal>* fcs, // source
      const FieldFace<Scal>* ffv, // volume flux
      double t, double dt, std::shared_ptr<Par> par)
      : ConvectionDiffusionScalar<M>(t, dt, fcr, ffd, fcs, ffv)
      , m(m) , mfc_(mfc) , mcc_(mcc) , par(par)
  {
    fc_field_.time_curr = fcu;
    fc_field_.time_prev = fc_field_.time_curr;
  }
  void StartStep() override {
    this->GetIter();
    if (IsNan(fc_field_.time_curr)) {
      throw std::runtime_error("NaN initial field");
    }
    fc_field_.iter_curr = fc_field_.time_curr;
    Scal ge = par->guessextra;
    for (auto idxcell : m.Cells()) {
      fc_field_.iter_curr[idxcell] +=
          (fc_field_.time_curr[idxcell] - fc_field_.time_prev[idxcell]) * ge;
    }
  }
  void MakeIteration() override {
    auto sem = m.GetSem("convdiff-iter");

    auto& fc_prev = fc_field_.iter_prev;
    auto& fc_curr = fc_field_.iter_curr;
    if (sem("assemble")) {
      fc_prev = fc_curr;

      fc_grad_ = Gradient(Interpolate(fc_prev, mfc_, m), m);

      InterpolationInnerFaceSecondUpwindDeferred<M, Expr>
      value_inner(m, *this->ffv_, fc_prev, fc_grad_);

      InterpolationBoundaryFaceNearestCell<M, Expr>
      value_boundary(m, mfc_);

      DerivativeInnerFacePlain<M, Expr>
      derivative_inner(m);

      DerivativeBoundaryFacePlain<M, Expr>
      derivative_boundary(m, mfc_);

      // Compute convective fluxes
			// all inner
      ff_cflux_.Reinit(m, Expr());
      for (IdxFace f : m.Faces()) {
        Expr e = value_inner.GetExpression(f);
        ff_cflux_[f] = e * (*this->ffv_)[f];
      }
			// overwrite with bc
      for (auto it = mfc_.cbegin(); it != mfc_.cend(); ++it) {
        IdxFace f = it->GetIdx();
        Expr e = value_boundary.GetExpression(f);
        ff_cflux_[f] = e * (*this->ffv_)[f];
      }

      // Compute diffusive fluxes
      // all inner
      ff_dflux_.Reinit(m, Expr());
      for (IdxFace f : m.Faces()) {
        Expr e = derivative_inner.GetExpression(f);
        ff_dflux_[f] = e *
            (-(*this->ffd_)[f]) * m.GetArea(f);
      }
			// overwrite with bc
      for (auto it = mfc_.cbegin(); it != mfc_.cend(); ++it) {
        IdxFace f = it->GetIdx();
				Expr e = derivative_boundary.GetExpression(f);
        ff_dflux_[f] = e *
            (-(*this->ffd_)[f]) * m.GetArea(f);
      }

      // Assemble the system
      fc_system_.Reinit(m);
      const Scal relax = par->relax;
      const bool second = par->second;
      for (IdxCell idxcell : m.Cells()) {
        Expr& eqn = fc_system_[idxcell];
        Expr cflux_sum;
        for (size_t i = 0; i < m.GetNumNeighbourFaces(idxcell); ++i) {
          IdxFace idxface = m.GetNeighbourFace(idxcell, i);
          cflux_sum += ff_cflux_[idxface] * m.GetOutwardFactor(idxcell, i);
        }

        Expr dflux_sum;
        for (size_t i = 0; i < m.GetNumNeighbourFaces(idxcell); ++i) {
          IdxFace idxface = m.GetNeighbourFace(idxcell, i);
          dflux_sum += ff_dflux_[idxface] * m.GetOutwardFactor(idxcell, i);
        }

        auto dt = this->GetTimeStep();
        auto coeffs = GetDerivativeApproxCoeffs(
            0., {-2. * dt, -dt, 0.}, second ? 0 : 1);

        Expr unsteady;
        unsteady.InsertTerm(coeffs[2], idxcell);
        unsteady.SetConstant(
            coeffs[0] * fc_field_.time_prev[idxcell] +
            coeffs[1] * fc_field_.time_curr[idxcell]);

        eqn = (unsteady + cflux_sum / m.GetVolume(idxcell)) *
              ((*this->fcr_)[idxcell]) +
              dflux_sum / m.GetVolume(idxcell) -
              Expr((*this->fcs_)[idxcell]);

        // Convert to delta-form
        eqn.SetConstant(eqn.Evaluate(fc_prev));

        // Apply under-relaxation
        eqn[eqn.Find(idxcell)].coeff /= relax;
      }

      // Fill halo cells with u=0 equations
      /*
      for (IdxCell idxcell : m.AllCells()) {
        Expr& eqn = fc_system_[idxcell];
        if (eqn.size() == 0) {
          eqn.Clear();
          eqn.InsertTerm(1., idxcell);
          eqn.SetConstant(0.);
        }
      }
      */

      // Include cell conditions for velocity
      for (auto it = mcc_.cbegin(); it != mcc_.cend(); ++it) {
        IdxCell c(it->GetIdx());
        ConditionCell* cond = it->GetValue().get();
        auto& eqn = fc_system_[c];
        if (auto cond_value = dynamic_cast<ConditionCellValue<Scal>*>(cond)) {
          eqn.Clear();
          // TODO: Revise dt coefficient for fixed-value cell condition
          eqn.InsertTerm(1. / this->GetTimeStep(), c);
          eqn.SetConstant((fc_prev[c] - cond_value->GetValue()) /
              this->GetTimeStep());
        }
      }
    }

    if (sem("solve")) {
      auto l = ConvertLs(fc_system_, lsa_, lsb_, lsx_, m);
      using T = typename M::LS::T;
      l.t = T::gen;
      m.Solve(l);
    }

    if (sem("applycomm")) {
      fc_corr_.Reinit(m);
      size_t i = 0;
      for (auto c : m.Cells()) {
        fc_corr_[c] = lsx_[i++];
      }
      assert(i == lsx_.size());
      for (auto idxcell : m.Cells()) {
        fc_curr[idxcell] = fc_prev[idxcell] + fc_corr_[idxcell];
      }
      m.Comm(&fc_curr);
      this->IncIter();
    }
  }
  void FinishStep() override {
    fc_field_.time_prev = fc_field_.time_curr;
    fc_field_.time_curr = fc_field_.iter_curr;
    if (IsNan(fc_field_.time_curr)) {
      throw std::runtime_error("NaN field");
    }
    this->IncTime();
  }
  double GetError() const override {
    if (this->GetIter() == 0) {
      return 1.;
    }
    return CalcDiff(fc_field_.iter_curr, fc_field_.iter_prev, m);
  }
  const FieldCell<Scal>& GetField() override {
    return fc_field_.time_curr;
  }
  const FieldCell<Scal>& GetField(Layers layer) override {
    return fc_field_.Get(layer);
  }
  void CorrectField(Layers layer,
                    const FieldCell<Scal>& fc_corr) override {
    auto sem = m.GetSem("convdiff-corr");
    if (sem("applycomm")) {
      auto& fc_field_layer = fc_field_.Get(layer);
      for (auto idxcell : m.Cells()) {
        fc_field_layer[idxcell] += fc_corr[idxcell];
      }
      m.Comm(&fc_field_layer);
    }
  }
  const FieldCell<Expr>& GetEquations() override {
    return fc_system_;
  }
};


} // namespace solver
