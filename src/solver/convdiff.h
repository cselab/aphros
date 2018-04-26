#pragma once

#include "solver.h"

namespace solver {

template <class Mesh>
class ConvectionDiffusionScalar : public UnsteadyIterativeSolver {
  using Scal = typename Mesh::Scal;
  static constexpr size_t dim = Mesh::dim;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

 protected:
  const FieldCell<Scal>* p_fc_scaling_;  // density
  const FieldFace<Scal>* p_ff_diffusion_rate_;  // dynamic viscosity
  const FieldCell<Scal>* p_fc_source_;
  const FieldFace<Scal>* p_ff_vol_flux_;

 public:
  ConvectionDiffusionScalar(
      double t, double dt,
      const FieldCell<Scal>* p_fc_scaling,
      const FieldFace<Scal>* p_ff_diffusion_rate,
      const FieldCell<Scal>* p_fc_source,
      const FieldFace<Scal>* p_ff_vol_flux)
      : UnsteadyIterativeSolver(t, dt)
      , p_fc_scaling_(p_fc_scaling)
      , p_ff_diffusion_rate_(p_ff_diffusion_rate)
      , p_fc_source_(p_fc_source)
      , p_ff_vol_flux_(p_ff_vol_flux)
  {

  }
  virtual const FieldCell<Scal>& GetField() = 0;
  virtual const FieldCell<Scal>& GetField(Layers layer) = 0;
  virtual void CorrectField(Layers layer,
                            const FieldCell<Scal>& fc_corr) = 0;
  virtual const FieldCell<Expr>& GetEquations() = 0;
};


} // namespace solver
