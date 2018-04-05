#pragma once

#include <exception>
#include <fstream>

#include "mesh.hpp"
#include "linear.hpp"
#include "solver.hpp"

namespace solver {

template <class M>
geom::FieldCell<typename M::Vect> GetDeformingVelocity(const M& m) {
  using Vect = typename M::Vect;
  geom::FieldCell<Vect> r(m, 0);
  for (auto c : m.Cells()) {
    auto x = m.GetCenter(c);
    r[c][0] = -std::cos(x[1]) * std::sin(x[0]);
    r[c][1] = std::cos(x[0]) * std::sin(x[1]);
  }
  return r;
}

template <class Mesh>
class AdvectionSolver : public UnsteadyIterativeSolver {
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;

 protected:
  const geom::FieldFace<Scal>* p_ff_velocity_;
  const geom::FieldCell<Scal>* p_fc_source_;

 public:
  AdvectionSolver(double t, double dt,
                  const geom::FieldFace<Scal>* p_ff_velocity,
                  const geom::FieldCell<Scal>* p_fc_source)
      : UnsteadyIterativeSolver(t, dt)
      , p_ff_velocity_(p_ff_velocity)
      , p_fc_source_(p_fc_source)
  {}
  virtual const geom::FieldCell<Scal>& GetField() = 0;
};

template <class Mesh>
class AdvectionSolverExplicit : public AdvectionSolver<Mesh> {
  Mesh& m;
  static constexpr size_t dim = Mesh::dim;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  LayersData<geom::FieldCell<Scal>> fc_u_;
  using IdxCell = geom::IdxCell;
  using IdxFace = geom::IdxFace;
  using IdxNode = geom::IdxNode;
  geom::MapFace<std::shared_ptr<ConditionFace>> mf_u_cond_;
  // Common buffers:
  geom::FieldFace<Vect> ff_velocity_;
  geom::FieldFace<Scal> ff_volume_flux_;
  geom::FieldFace<Scal> ff_flux_;
  geom::FieldFace<Scal> ff_u_;
  geom::FieldFace<Scal> sharp_af_;
  geom::FieldCell<Vect> sharp_gc_;
  geom::FieldFace<Vect> sharp_gf_;
  Scal sharp_;
  Scal sharpo_;
  Scal sharp_max_;

 public:
  AdvectionSolverExplicit(
      Mesh& m,
      const geom::FieldCell<Scal>& fc_u_initial,
      const geom::MapFace<std::shared_ptr<ConditionFace>>& mf_u_cond_,
      const geom::FieldFace<Scal>* p_fn_velocity,
      const geom::FieldCell<Scal>* p_fc_source,
      double time, double time_step,
      Scal sharp, Scal sharpo, Scal sharp_max)
      : AdvectionSolver<Mesh>(
          time, time_step, p_fn_velocity, p_fc_source)
      , m(m)
      , mf_u_cond_(mf_u_cond_)
      , ff_volume_flux_(m)
      , ff_flux_(m)
      , sharp_(sharp)
      , sharpo_(sharpo)
      , sharp_max_(sharp_max)
  {
    fc_u_.time_curr = fc_u_initial;
  }
  void StartStep() override {
    this->ClearIter();
    fc_u_.time_prev = fc_u_.time_curr;
    fc_u_.iter_curr = fc_u_.time_prev;
  }
  void MakeIteration() override {
    auto sem = m.GetSem();
    if (sem()) {
      auto& prev = fc_u_.iter_prev;
      auto& curr = fc_u_.iter_curr;
      prev = curr;

      ff_volume_flux_ = *(this->p_ff_velocity_);

      ff_u_ = InterpolateSuperbee(
          prev,
          Gradient(Interpolate(prev, mf_u_cond_, m), m),
          mf_u_cond_, ff_volume_flux_, m);

      for (auto idxface : m.Faces()) {
        ff_flux_[idxface] = ff_u_[idxface] * ff_volume_flux_[idxface];
      }

      // Interface sharpening
      // append to ff_flux_
      if (std::abs(sharp_) != 0.) {
        // zero-derivative bc for Vect
        geom::MapFace<std::shared_ptr<ConditionFace>> mfvz;
        for (auto it : mf_u_cond_) {
          IdxFace i = it.GetIdx();
          mfvz[i] = std::make_shared<
              ConditionFaceDerivativeFixed<Vect>>(
                  Vect(0), it.GetValue()->GetNci());
        }
        auto& af = sharp_af_;
        auto& gc = sharp_gc_;
        auto& gf = sharp_gf_;
        af = Interpolate(curr, mf_u_cond_, m);
        gc = Gradient(af, m);
        gf = Interpolate(gc, mfvz, m);
        for (auto i : m.Faces()) {
          // normal to interface
          const Vect ni = gf[i] / (gf[i].norm() + 1e-6); 
          const Vect n = m.GetNormal(i);
          const Vect s = m.GetSurface(i);
          const Scal a = m.GetArea(i);
          //const Scal nf = n.dot(m.GetNormal(i));
          //const Scal uf = ff_volume_flux_[i];
          //const Scal uf = 1.;
          //const Scal am = sharp_max_;
          //const Scal eh = sharp_ * m.GetArea(i);
          //ff_flux_[i] -= sharpo_ * std::abs(uf * nf) * nf * 
          //    (eh * gf[i].norm() - af[i] * (1. - af[i] / am));
          ff_flux_[i] += 
            sharp_ * af[i] * (1. - af[i]) * ni.dot(s)
            - sharpo_ * gf[i].dot(s);
        }
      }

      for (auto idxcell : m.Cells()) {
        Scal flux_sum = 0.;
        for (size_t i = 0; i < m.GetNumNeighbourFaces(idxcell); ++i) {
          IdxFace idxface = m.GetNeighbourFace(idxcell, i);
          flux_sum += ff_flux_[idxface] * m.GetOutwardFactor(idxcell, i);
        }

        curr[idxcell] = fc_u_.time_prev[idxcell] -
            this->GetTimeStep() / m.GetVolume(idxcell) * flux_sum +
            this->GetTimeStep() * (*this->p_fc_source_)[idxcell];
      }
      m.Comm(&curr);

      this->IncIter();
    }
  }
  void FinishStep() override {
    fc_u_.time_curr = fc_u_.iter_curr;
    this->IncTime();
  }
  double GetError() const override {
    if (this->GetError() == 0) {
      return 1.;
    }
    return CalcDiff<geom::FieldCell<Scal>, Mesh>(
        fc_u_.iter_curr, fc_u_.iter_prev, m);
  }
  const geom::FieldCell<Scal>& GetField() override {
    return fc_u_.time_curr;
  }
};

} // namespace solver
