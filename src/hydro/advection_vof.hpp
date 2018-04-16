#pragma once

#include <exception>
#include <fstream>

#include "mesh.hpp"
#include "linear.hpp"
#include "solver.hpp"
#include "advection.hpp"

namespace solver {

using namespace geom;

template <class M>
class Vof : public AdvectionSolver<M> {
  using Mesh = M;
  using P = AdvectionSolver<Mesh>;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  static constexpr size_t dim = Mesh::dim;

  using P::m;
  using P::ffv_;
  LayersData<FieldCell<Scal>> fc_u_;
  MapFace<std::shared_ptr<ConditionFace>> mf_u_cond_;

  // Equation of reconstructed plane would be
  // ((x-xc)/h).dot(n) = alpha
  // where xc -- cell center, h - cell size, n -- unit normal to plane
  FieldCell<Scal> fc_a_; // alpha (plane constant)
  FieldCell<Vect> fc_n_; // n (normal to plane)

 public:
  struct Par {
  };
  Par* par;
  Vof(
      Mesh& m,
      const FieldCell<Scal>& fc_u_initial,
      const MapFace<std::shared_ptr<ConditionFace>>& mf_u_cond_,
      const FieldFace<Scal>* ff_volume_flux,
      const FieldCell<Scal>* fc_source,
      double t, double dt, Par* par)
      : AdvectionSolver<Mesh>(t, dt, m, ff_volume_flux, fc_source)
      , mf_u_cond_(mf_u_cond_)
      , par(par)
      , fc_a_(m, 0)
      , fc_n_(m, Vect(0))
  {
    fc_u_.time_curr = fc_u_initial;
  }
  void StartStep() override {
    this->ClearIter();
    fc_u_.time_prev = fc_u_.time_curr;
    fc_u_.iter_curr = fc_u_.time_prev;
  }
  void MakeIteration() override {
    auto sem = m.GetSem("iter");
    if (sem("init")) {
      auto& curr = fc_u_.iter_curr;
      const Scal dt = this->GetTimeStep();
      for (auto c : m.Cells()) {
        curr[c] = fc_u_.time_prev[c] +  // previous time step
            dt * (*this->fcs_)[c]; // source
      }
    }
    const bool split = true;
    for (size_t d = 0; d < (split ? dim : 1); ++d) {
      if (sem("adv")) {
        auto& curr = fc_u_.iter_curr;

        

        auto& ffv = *ffv_; // [f]ield [f]ace [v]olume flux
        m.Comm(&curr);
      }
    }

    if (sem("stat")) {
      this->IncIter();
    }
  }
  void FinishStep() override {
    fc_u_.time_curr = fc_u_.iter_curr;
    this->IncTime();
  }
  const FieldCell<Scal>& GetField(Layers l) override {
    return fc_u_.Get(l);
  }
  const FieldCell<Scal>& GetAlpha() {
    return fc_a_;
  }
  const FieldCell<Vect>& GetNormal() {
    return fc_n_;
  }
  using P::GetField;
};

} // namespace solver
