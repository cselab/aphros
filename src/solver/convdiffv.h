#pragma once

#include "convdiff.h"

namespace solver {

template <class Mesh>
class ConvectionDiffusion : public UnsteadyIterativeSolver {
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  static constexpr size_t dim = Mesh::dim;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

 protected:
  FieldCell<Scal>* p_fc_density_;
  // TODO: rename to dynamic viscosity
  FieldFace<Scal>* p_ff_kinematic_viscosity_; 
  FieldCell<Vect>* p_fc_force_;
  FieldFace<Scal>* p_ff_vol_flux_;

 public:
  ConvectionDiffusion(double time, double time_step,
                      FieldCell<Scal>* p_fc_density,
                      FieldFace<Scal>* p_ff_kinematic_viscosity,
                      FieldCell<Vect>* p_fc_force,
                      FieldFace<Scal>* p_ff_vol_flux
                      )
      : UnsteadyIterativeSolver(time, time_step)
      , p_fc_density_(p_fc_density)
      , p_ff_kinematic_viscosity_(p_ff_kinematic_viscosity)
      , p_fc_force_(p_fc_force)
      , p_ff_vol_flux_(p_ff_vol_flux)
  {

  }
  virtual const FieldCell<Vect>& GetVelocity() = 0;
  virtual const FieldCell<Vect>& GetVelocity(Layers layer) = 0;
  virtual void CorrectVelocity(Layers layer,
                               const FieldCell<Vect>& fc_corr) = 0;
  virtual const FieldCell<Expr>& GetVelocityEquations(size_t comp) = 0;
};


} // namespace solver
