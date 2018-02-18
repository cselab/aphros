/*
 * conv_diff.hpp
 *
 *  Created on: May 5, 2016
 *      Author: Petr Karnakov (petr.karnakov@gmail.com)
 */

#pragma once

#include "mesh.hpp"
#include "linear.hpp"
#include <exception>
#include "solver.hpp"
#include "../control/metrics.hpp"
#include <cmath>
#include <sstream>

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
      double time, double time_step,
      const geom::FieldCell<Scal>* p_fc_scaling,
      const geom::FieldFace<Scal>* p_ff_diffusion_rate,
      const geom::FieldCell<Scal>* p_fc_source,
      const geom::FieldFace<Scal>* p_ff_vol_flux,
      double convergence_tolerance,
      size_t num_iterations_limit)
      : UnsteadyIterativeSolver(time, time_step, convergence_tolerance,
                                num_iterations_limit)
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
  const Mesh& mesh;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  static constexpr size_t dim = Mesh::dim;
  LayersData<geom::FieldCell<Scal>> fc_field_;
  using IdxCell = geom::IdxCell;
  using IdxFace = geom::IdxFace;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

  geom::MapFace<std::shared_ptr<ConditionFace>> mf_cond_;
  geom::MapCell<std::shared_ptr<ConditionCell>> mc_cond_;
  std::shared_ptr<LinearSolver<Scal, IdxCell, Expr>> linear_;
  Scal relaxation_factor_;
  bool time_second_order_;
  Scal guess_extrapolation_;

  // Common buffers:
  geom::FieldFace<Expr> ff_cflux_;
  geom::FieldFace<Expr> ff_dflux_;
  geom::FieldCell<Expr> fc_system_;
  geom::FieldCell<Scal> fc_corr_;
  geom::FieldCell<Vect> fc_grad_;

 public:
  ConvectionDiffusionScalarImplicit(
      const Mesh& mesh,
      const geom::FieldCell<Scal>& fc_initial,
      const geom::MapFace<std::shared_ptr<ConditionFace>>&
      mf_cond,
      const geom::MapCell<std::shared_ptr<ConditionCell>>&
      mc_cond,
      Scal relaxation_factor,
      const geom::FieldCell<Scal>* p_fc_scaling,
      const geom::FieldFace<Scal>* p_ff_diffusion_rate,
      const geom::FieldCell<Scal>* p_fc_source,
      const geom::FieldFace<Scal>* p_ff_vol_flux,
      double time, double time_step,
      const LinearSolverFactory& linear_factory,
      double convergence_tolerance,
      size_t num_iterations_limit,
      bool time_second_order = true,
      Scal guess_extrapolation = 0.)
      : ConvectionDiffusionScalar<Mesh>(
          time, time_step, p_fc_scaling, p_ff_diffusion_rate, p_fc_source,
          p_ff_vol_flux,
          convergence_tolerance, num_iterations_limit)
      , mesh(mesh)
      , mf_cond_(mf_cond)
      , mc_cond_(mc_cond)
      , relaxation_factor_(relaxation_factor)
      , time_second_order_(time_second_order)
      , guess_extrapolation_(guess_extrapolation)
  {
    fc_field_.time_curr = fc_initial;
    fc_field_.time_prev = fc_field_.time_curr;

    linear_ = linear_factory.Create<Scal, IdxCell, Expr>();
  }
  void StartStep() override {
    this->ClearIterationCount();
    if (IsNan(fc_field_.time_curr)) {
      throw std::string("NaN initial field");
    }
    fc_field_.iter_curr = fc_field_.time_curr;
    for (auto idxcell : mesh.Cells()) {
      fc_field_.iter_curr[idxcell] +=
          (fc_field_.time_curr[idxcell] - fc_field_.time_prev[idxcell]) *
          guess_extrapolation_;
    }
  }
  void MakeIteration() override {
    auto& fc_prev = fc_field_.iter_prev;
    auto& fc_curr = fc_field_.iter_curr;
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
    ff_cflux_.Reinit(mesh, Expr());
#pragma omp parallel for
    for (IntIdx i = 0; i < static_cast<IntIdx>(mesh.Faces().size()); ++i) {
      IdxFace idxface(i);
      if (!mesh.IsExcluded(idxface)) {
        Expr value_expr, derivative_expr;
        if (mesh.IsInner(idxface)) {
          value_expr = value_inner.GetExpression(idxface);
        } else {
          value_expr = value_boundary.GetExpression(idxface);
        }
        ff_cflux_[idxface] =
            value_expr * (*this->p_ff_vol_flux_)[idxface];
      }
    }

    // Compute diffusive fluxes
    ff_dflux_.Reinit(mesh, Expr());
#pragma omp parallel for
    for (IntIdx i = 0; i < static_cast<IntIdx>(mesh.Faces().size()); ++i) {
      IdxFace idxface(i);
      if (!mesh.IsExcluded(idxface)) {
        Expr value_expr, derivative_expr;
        if (mesh.IsInner(idxface)) {
          derivative_expr = derivative_inner.GetExpression(idxface);
        } else {
          derivative_expr = derivative_boundary.GetExpression(idxface);
        }
        ff_dflux_[idxface] = derivative_expr *
            (-(*this->p_ff_diffusion_rate_)[idxface]) * mesh.GetArea(idxface);
      }
    }

    // Assemble the system
    fc_system_.Reinit(mesh);
#pragma omp parallel for
    for (IntIdx i = 0; i < static_cast<IntIdx>(mesh.Cells().size()); ++i) {
      IdxCell idxcell(i);
      Expr& eqn = fc_system_[idxcell];
      if (!mesh.IsExcluded(idxcell)) {
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
            0., {-2. * dt, -dt, 0.}, time_second_order_ ? 0 : 1);

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
        eqn[eqn.Find(idxcell)].coeff /= relaxation_factor_;
      } else {
        eqn.Clear();
        eqn.InsertTerm(1., idxcell);
        eqn.SetConstant(0.);
      }
    }

    // Account for cell conditions for velocity
    for (auto it = mc_cond_.cbegin();
        it != mc_cond_.cend(); ++it) {
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

    fc_corr_ = linear_->Solve(fc_system_);
    for (auto idxcell : mesh.Cells()) {
      fc_curr[idxcell] = fc_prev[idxcell] + fc_corr_[idxcell];
    }

    this->IncIterationCount();
  }
  void FinishStep() override {
    fc_field_.time_prev = fc_field_.time_curr;
    fc_field_.time_curr = fc_field_.iter_curr;
    if (IsNan(fc_field_.time_curr)) {
      throw std::string("NaN field");
    }
    this->IncTime();
  }
  double GetConvergenceIndicator() const override {
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
    auto& fc_field_layer = fc_field_.Get(layer);
    for (auto idxcell : mesh.Cells()) {
      fc_field_layer[idxcell] += fc_corr[idxcell];
    }
  }
  const geom::FieldCell<Expr>& GetEquations() override {
    return fc_system_;
  }
};


} // namespace solver
