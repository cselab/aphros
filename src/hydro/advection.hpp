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
  using P::ffv_;
  LayersData<FieldCell<Scal>> fc_u_;
  MapFace<std::shared_ptr<ConditionFace>> mf_u_cond_;
  // Common buffers:
  FieldFace<Vect> ff_velocity_;
  FieldFace<Scal> ff_volume_flux_;
  FieldFace<Scal> ff_flux_;
  FieldFace<Scal> ff_u_;
  FieldCell<Scal> sharp_as_;
  FieldFace<Scal> sharp_af_;
  FieldCell<Vect> sharp_gc_;
  FieldFace<Vect> sharp_gf_;
  Scal sharp_;
  Scal sharpo_;
  Scal sharp_max_;

 public:
  struct Par {
    Scal sharp = 0.;
    Scal sharpo = 0.;
    Scal sharp_max = 1.;
    bool split = false;
  };
  Par* par;
  AdvectionSolverExplicit(
      Mesh& m,
      const FieldCell<Scal>& fc_u_initial,
      const MapFace<std::shared_ptr<ConditionFace>>& mf_u_cond_,
      const FieldFace<Scal>* p_fn_velocity,
      const FieldCell<Scal>* p_fc_source,
      double t, double dt, Par* par)
      : AdvectionSolver<Mesh>(t, dt, m, p_fn_velocity, p_fc_source)
      , mf_u_cond_(mf_u_cond_)
      , ff_volume_flux_(m)
      , ff_flux_(m)
      , par(par)
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
    for (size_t d = 0; d < (par->split ? dim : 1); ++d) {
      if (sem("adv")) {
        auto& curr = fc_u_.iter_curr;

        auto& ffv = *ffv_; // [f]ield [f]ace [v]olume flux
        
        ff_u_ = InterpolateSuperbee(
            curr,
            Gradient(Interpolate(curr, mf_u_cond_, m), m),
            mf_u_cond_, ffv, m);

        for (auto f : m.Faces()) {
          ff_flux_[f] = ff_u_[f] * ffv[f];
        }

        if (par->split) {
          Vect vd(0);
          vd[d] = 1.;

          // Weighten by n.dot(vd)
          for (auto f : m.Faces()) {
            auto n = m.GetNormal(f);
            ff_flux_[f] *= (n.dot(vd)) / n.norm1();
          }
        }

        const Scal dt = this->GetTimeStep();
        for (auto c : m.Cells()) {
          Scal s = 0.; // sum of fluxes
          for (auto q : m.Nci(c)) {
            IdxFace f = m.GetNeighbourFace(c, q);
            s += ff_flux_[f] * m.GetOutwardFactor(c, q);
          }

          curr[c] += dt / m.GetVolume(c) * (-s); // fluxes
        }
        m.Comm(&curr);

        sharp_as_ = fc_u_.iter_curr;
      }
    }

    //if (sem.Nested("smooth")) {
      //Smoothen(sharp_as_, mf_u_cond_, m, 5);
    //}

    // Interface sharpening
    if (std::abs(par->sharp) != 0. && sem("sharp")) {
      auto& ac = fc_u_.iter_curr;

      // zero-derivative bc for Vect
      MapFace<std::shared_ptr<ConditionFace>> mfvz;
      for (auto it : mf_u_cond_) {
        IdxFace f = it.GetIdx();
        mfvz[f] = std::make_shared<
            ConditionFaceDerivativeFixed<Vect>>(
                Vect(0), it.GetValue()->GetNci());
      }

      // smooth
      //auto asf = Interpolate(sharp_as_, mf_u_cond_, m);
      //auto gsc = Gradient(asf, m);
      //auto gsf = Interpolate(gsc, mfvz, m);

      auto& af = sharp_af_;
      auto& gc = sharp_gc_;
      auto& gf = sharp_gf_;
      af = Interpolate(ac, mf_u_cond_, m);
      gc = Gradient(af, m);
      gf = Interpolate(gc, mfvz, m);
      const Scal kc = par->sharp;
      const Scal kd = par->sharp * par->sharpo;

      for (auto f : m.Faces()) {
        const Vect g = gf[f];  // gradient on face
        //const Vect gs = gsf[f];  // gradient on face
        const Vect n = g / (g.norm() + 1e-6);  // normal to interface
        //const Vect ns = gs / (gs.norm() + 1e-6); 
        const Vect s = m.GetSurface(f); // surface vector to face
        IdxCell cm = m.GetNeighbourCell(f, 0);
        IdxCell cp = m.GetNeighbourCell(f, 1);
        Scal am = ac[cm];
        Scal ap = ac[cp];
        ap = std::max(0., std::min(1., ap));
        am = std::max(0., std::min(1., am));
        ff_flux_[f] = 
          kc * (ap > am 
             ?  (1. - (ap - am)) * (1. - ap) * am
             :  -(1. - (am - ap)) * (1. - am) * ap)
          -kd * (ap - am);
        ff_flux_[f] *= std::abs(n.dot(s));
      }

      for (auto c : m.Cells()) {
        Scal s = 0.; // sum of fluxes
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetNeighbourFace(c, q);
          s += ff_flux_[f] * m.GetOutwardFactor(c, q);
        }

        const Scal dt = this->GetTimeStep();
        ac[c] += dt / m.GetVolume(c) * (-s);
      }
      m.Comm(&ac);

    }
    if (sem("stat")) {
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
