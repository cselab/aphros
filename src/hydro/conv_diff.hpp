/*
 * conv_diff.hpp
 *
 *  Created on: May 5, 2016
 *      Author: Petr Karnakov (petr.karnakov@gmail.com)
 */

#pragma once

#include <cmath>
#include <sstream>
#include <stdexcept>

#include "mesh.hpp"
#include "linear.hpp"
#include "solver.hpp"

namespace solver {

template <class Mesh>
class ConvectionDiffusionScalar : public UnsteadyIterativeSolver {
  using Scal = typename Mesh::Scal;
  static constexpr size_t dim = Mesh::dim;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

 protected:
  const geom::FieldCell<Scal>* p_fc_scaling_;  // density
  const geom::FieldFace<Scal>* p_ff_diffusion_rate_;  // dynamic viscosity
  const geom::FieldCell<Scal>* p_fc_source_;
  const geom::FieldFace<Scal>* p_ff_vol_flux_;

 public:
  ConvectionDiffusionScalar(
      double t, double dt,
      const geom::FieldCell<Scal>* p_fc_scaling,
      const geom::FieldFace<Scal>* p_ff_diffusion_rate,
      const geom::FieldCell<Scal>* p_fc_source,
      const geom::FieldFace<Scal>* p_ff_vol_flux)
      : UnsteadyIterativeSolver(t, dt)
      , p_fc_scaling_(p_fc_scaling)
      , p_ff_diffusion_rate_(p_ff_diffusion_rate)
      , p_fc_source_(p_fc_source)
      , p_ff_vol_flux_(p_ff_vol_flux)
  {

  }
  virtual const geom::FieldCell<Scal>& GetField() = 0;
  virtual const geom::FieldCell<Scal>& GetField(Layers layer) = 0;
  virtual void CorrectField(Layers layer,
                            const geom::FieldCell<Scal>& fc_corr) = 0;
  virtual const geom::FieldCell<Expr>& GetEquations() = 0;
};

template <class Mesh>
class ConvectionDiffusionScalarImplicit :
    public ConvectionDiffusionScalar<Mesh> {
  Mesh& mesh;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  static constexpr size_t dim = Mesh::dim;
  LayersData<geom::FieldCell<Scal>> fc_field_;
  using IdxCell = geom::IdxCell;
  using IdxFace = geom::IdxFace;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

  geom::MapFace<std::shared_ptr<ConditionFace>> mf_cond_;
  geom::MapCell<std::shared_ptr<ConditionCell>> mc_cond_;

  // Common buffers:
  geom::FieldFace<Expr> ff_cflux_;
  geom::FieldFace<Expr> ff_dflux_;
  geom::FieldCell<Expr> fc_system_;
  geom::FieldCell<Scal> fc_corr_;
  geom::FieldCell<Vect> fc_grad_;
  // LS
  std::vector<Scal> lsa_, lsb_, lsx_;

 public:
  struct Par {
    Scal relax = 1.;      // relaxation factor [0,1] (1 -- no relaxation)
    Scal guessextra = 0.; // next iteration guess extrapolation weight [0,1]
    bool second = true; // second order in time
  };
  Par* par;
  ConvectionDiffusionScalarImplicit(
      Mesh& mesh,
      const geom::FieldCell<Scal>& fc_initial,
      const geom::MapFace<std::shared_ptr<ConditionFace>>& mf_cond,
      const geom::MapCell<std::shared_ptr<ConditionCell>>& mc_cond,
      const geom::FieldCell<Scal>* p_fc_scaling,
      const geom::FieldFace<Scal>* p_ff_diffusion_rate,
      const geom::FieldCell<Scal>* p_fc_source,
      const geom::FieldFace<Scal>* p_ff_vol_flux,
      double t, double dt, Par* par)
      : ConvectionDiffusionScalar<Mesh>(
          t, dt, 
          p_fc_scaling, p_ff_diffusion_rate, p_fc_source, p_ff_vol_flux)
      , mesh(mesh)
      , mf_cond_(mf_cond)
      , mc_cond_(mc_cond)
      , par(par)
  {
    fc_field_.time_curr = fc_initial;
    fc_field_.time_prev = fc_field_.time_curr;
  }
  void StartStep() override {
    this->ClearIterationCount();
    if (IsNan(fc_field_.time_curr)) {
      throw std::runtime_error("NaN initial field");
    }
    fc_field_.iter_curr = fc_field_.time_curr;
    Scal ge = par->guessextra;
    for (auto idxcell : mesh.Cells()) {
      fc_field_.iter_curr[idxcell] +=
          (fc_field_.time_curr[idxcell] - fc_field_.time_prev[idxcell]) * ge;
    }
  }
  void MakeIteration() override {
    auto sem = mesh.GetSem("convdiff-iter");
    auto& m = mesh;

    auto& fc_prev = fc_field_.iter_prev;
    auto& fc_curr = fc_field_.iter_curr;
    if (sem("assemble")) {
      fc_prev = fc_curr;

      fc_grad_ = Gradient(Interpolate(fc_prev, mf_cond_, mesh), mesh);

      InterpolationInnerFaceSecondUpwindDeferred<Mesh, Expr>
      value_inner(mesh, *this->p_ff_vol_flux_, fc_prev, fc_grad_);

      InterpolationBoundaryFaceNearestCell<Mesh, Expr>
      value_boundary(mesh, mf_cond_);

      DerivativeInnerFacePlain<Mesh, Expr>
      derivative_inner(mesh);

      DerivativeBoundaryFacePlain<Mesh, Expr>
      derivative_boundary(mesh, mf_cond_);

      // Compute convective fluxes
			// all inner
      ff_cflux_.Reinit(mesh, Expr());
      for (IdxFace f : mesh.Faces()) {
        Expr e = value_inner.GetExpression(f);
        ff_cflux_[f] = e * (*this->p_ff_vol_flux_)[f];
      }
			// overwrite with bc
      for (auto it = mf_cond_.cbegin(); it != mf_cond_.cend(); ++it) {
        IdxFace f = it->GetIdx();
        Expr e = value_boundary.GetExpression(f);
        ff_cflux_[f] = e * (*this->p_ff_vol_flux_)[f];
      }

      // Compute diffusive fluxes
      // all inner
      ff_dflux_.Reinit(mesh, Expr());
      for (IdxFace f : mesh.Faces()) {
        Expr e = derivative_inner.GetExpression(f);
        ff_dflux_[f] = e *
            (-(*this->p_ff_diffusion_rate_)[f]) * mesh.GetArea(f);
      }
			// overwrite with bc
      for (auto it = mf_cond_.cbegin(); it != mf_cond_.cend(); ++it) {
        IdxFace f = it->GetIdx();
				Expr e = derivative_boundary.GetExpression(f);
        ff_dflux_[f] = e *
            (-(*this->p_ff_diffusion_rate_)[f]) * mesh.GetArea(f);
      }

      // Assemble the system
      fc_system_.Reinit(mesh);
      const Scal relax = par->relax;
      const bool second = par->second;
      for (IdxCell idxcell : mesh.Cells()) {
        Expr& eqn = fc_system_[idxcell];
        Expr cflux_sum;
        for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
          IdxFace idxface = mesh.GetNeighbourFace(idxcell, i);
          cflux_sum += ff_cflux_[idxface] * mesh.GetOutwardFactor(idxcell, i);
        }

        Expr dflux_sum;
        for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
          IdxFace idxface = mesh.GetNeighbourFace(idxcell, i);
          dflux_sum += ff_dflux_[idxface] * mesh.GetOutwardFactor(idxcell, i);
        }

        auto dt = this->GetTimeStep();
        auto coeffs = GetDerivativeApproxCoeffs(
            0., {-2. * dt, -dt, 0.}, second ? 0 : 1);

        Expr unsteady;
        unsteady.InsertTerm(coeffs[2], idxcell);
        unsteady.SetConstant(
            coeffs[0] * fc_field_.time_prev[idxcell] +
            coeffs[1] * fc_field_.time_curr[idxcell]);

        eqn = (unsteady + cflux_sum / mesh.GetVolume(idxcell)) *
              ((*this->p_fc_scaling_)[idxcell]) +
              dflux_sum / mesh.GetVolume(idxcell) -
              Expr((*this->p_fc_source_)[idxcell]);

        // Convert to delta-form
        eqn.SetConstant(eqn.Evaluate(fc_prev));

        // Apply under-relaxation
        eqn[eqn.Find(idxcell)].coeff /= relax;
      }

      // Fill halo cells with u=0 equations
      /*
      for (IdxCell idxcell : mesh.AllCells()) {
        Expr& eqn = fc_system_[idxcell];
        if (eqn.size() == 0) {
          eqn.Clear();
          eqn.InsertTerm(1., idxcell);
          eqn.SetConstant(0.);
        }
      }
      */

      // Include cell conditions for velocity
      for (auto it = mc_cond_.cbegin(); it != mc_cond_.cend(); ++it) {
        IdxCell idxcell(it->GetIdx());
        ConditionCell* cond = it->GetValue().get();
        auto& eqn = fc_system_[idxcell];
        if (auto cond_value = dynamic_cast<ConditionCellValue<Scal>*>(cond)) {
          eqn.Clear();
          // TODO: Revise dt coefficient for fixed-value cell condition
          eqn.InsertTerm(1. / this->GetTimeStep(), idxcell);
          eqn.SetConstant(
              (fc_prev[idxcell] - cond_value->GetValue()) /
              this->GetTimeStep());
        }
      }
    }

    if (sem("solve")) {
      auto l = ConvertLs(fc_system_, lsa_, lsb_, lsx_, mesh);
      using T = typename Mesh::LS::T;
      l.t = T::gen;
      m.Solve(l);
    }

    if (sem("applycomm")) {
      fc_corr_.Reinit(mesh);
      size_t i = 0;
      for (auto c : m.Cells()) {
        fc_corr_[c] = lsx_[i++];
      }
      assert(i == lsx_.size());
      for (auto idxcell : mesh.Cells()) {
        fc_curr[idxcell] = fc_prev[idxcell] + fc_corr_[idxcell];
      }
      m.Comm(&fc_curr);
      this->IncIterationCount();
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
    if (this->GetIterationCount() == 0) {
      return 1.;
    }
    return CalcDiff(fc_field_.iter_curr, fc_field_.iter_prev, mesh);
  }
  const geom::FieldCell<Scal>& GetField() override {
    return fc_field_.time_curr;
  }
  const geom::FieldCell<Scal>& GetField(Layers layer) override {
    return fc_field_.Get(layer);
  }
  void CorrectField(Layers layer,
                    const geom::FieldCell<Scal>& fc_corr) override {
    auto sem = mesh.GetSem("convdiff-corr");
    if (sem("applycomm")) {
      auto& fc_field_layer = fc_field_.Get(layer);
      for (auto idxcell : mesh.Cells()) {
        fc_field_layer[idxcell] += fc_corr[idxcell];
      }
      mesh.Comm(&fc_field_layer);
    }
  }
  const geom::FieldCell<Expr>& GetEquations() override {
    return fc_system_;
  }
};


} // namespace solver
