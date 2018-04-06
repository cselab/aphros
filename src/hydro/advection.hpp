#pragma once

#include <exception>
#include <fstream>

#include "mesh.hpp"
#include "linear.hpp"
#include "solver.hpp"

namespace solver {

using namespace geom;

template <class M>
FieldCell<typename M::Vect> GetDeformingVelocity(const M& m) {
  using Vect = typename M::Vect;
  FieldCell<Vect> r(m, 0);
  for (auto c : m.Cells()) {
    auto x = m.GetCenter(c);
    r[c][0] = -std::cos(x[1]) * std::sin(x[0]);
    r[c][1] = std::cos(x[0]) * std::sin(x[1]);
  }
  return r;
}

template <class M>
class AdvectionSolver : public UnsteadyIterativeSolver {
  using Mesh = M;
  using Scal = typename Mesh::Scal;

 protected:
  Mesh& m;
  const FieldFace<Scal>* ffv_; // volume flux [velocity*area]
  const FieldCell<Scal>* fcs_; // source [value/time]

 public:
  AdvectionSolver(double t, double dt,
                  Mesh& m,
                  const FieldFace<Scal>* ffv /*volume flux*/,
                  const FieldCell<Scal>* fcs /*source*/)
      : UnsteadyIterativeSolver(t, dt)
      , m(m)
      , ffv_(ffv)
      , fcs_(fcs) {}
  virtual const FieldCell<Scal>& GetField(Layers) = 0;
  virtual const FieldCell<Scal>& GetField() {
    return GetField(Layers::time_curr);
  }
};

template <class M>
class AdvectionSolverExplicit : public AdvectionSolver<M> {
  using Mesh = M;
  using P = AdvectionSolver<Mesh>;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  static constexpr size_t dim = Mesh::dim;

  using P::m;
  LayersData<FieldCell<Scal>> fc_u_;
  MapFace<std::shared_ptr<ConditionFace>> mf_u_cond_;
  // Common buffers:
  FieldFace<Vect> ff_velocity_;
  FieldFace<Scal> ff_volume_flux_;
  FieldFace<Scal> ff_flux_;
  FieldFace<Scal> ff_u_;
  FieldFace<Scal> sharp_af_;
  FieldCell<Vect> sharp_gc_;
  FieldFace<Vect> sharp_gf_;
  Scal sharp_;
  Scal sharpo_;
  Scal sharp_max_;

 public:
  AdvectionSolverExplicit(
      Mesh& m,
      const FieldCell<Scal>& fc_u_initial,
      const MapFace<std::shared_ptr<ConditionFace>>& mf_u_cond_,
      const FieldFace<Scal>* p_fn_velocity,
      const FieldCell<Scal>* p_fc_source,
      double t, double dt,
      Scal sharp, Scal sharpo, Scal sharp_max)
      : AdvectionSolver<Mesh>(t, dt, m, p_fn_velocity, p_fc_source)
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

      ff_volume_flux_ = *(this->ffv_);

      ff_u_ = InterpolateSuperbee(
          prev,
          Gradient(Interpolate(prev, mf_u_cond_, m), m),
          mf_u_cond_, ff_volume_flux_, m);

      for (auto f : m.Faces()) {
        ff_flux_[f] = ff_u_[f] * ff_volume_flux_[f];
      }

      // Interface sharpening
      // append to ff_flux_
      if (std::abs(sharp_) != 0.) {
        // zero-derivative bc for Vect
        MapFace<std::shared_ptr<ConditionFace>> mfvz;
        for (auto it : mf_u_cond_) {
          IdxFace f = it.GetIdx();
          mfvz[f] = std::make_shared<
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

      for (auto c : m.Cells()) {
        Scal flux_sum = 0.;
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetNeighbourFace(c, q);
          flux_sum += ff_flux_[f] * m.GetOutwardFactor(c, q);
        }

        const Scal dt = this->GetTimeStep();
        curr[c] = fc_u_.time_prev[c] +
            dt / m.GetVolume(c) * (-flux_sum) + // fluxes
            dt * (*this->fcs_)[c]; // source
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
    return CalcDiff<FieldCell<Scal>, Mesh>(
        fc_u_.iter_curr, fc_u_.iter_prev, m);
  }
  const FieldCell<Scal>& GetField(Layers l) override {
    return fc_u_.Get(l);
  }
  using P::GetField;
};

} // namespace solver
