#pragma once

#include <exception>
#include <fstream>

#include "mesh.hpp"
#include "linear.hpp"
#include "solver.hpp"
#include "advection.hpp"

namespace solver {

using namespace geom;

// Volume fraction to line.
// Returns:
// a: line constant
// Equation of reconstructed plane would be
// ((x-xc)/h).dot(n) = a
// where 
// xc: cell center
// h: cell size
template <class Scal>
inline Scal GetLineA(const GVect<Scal, 3>& n /*unit normal*/, 
                     Scal u /*volume fraction*/) {
  using Vect = GVect<Scal, 3>;
  Scal alpha, mx, m2, v1;

  Scal nx = std::abs(n[0]); 
  Scal ny = std::abs(n[1]);
  if (nx > ny) {
    std::swap(nx, ny);
  }
  // nx < ny
  
  const Scal w = nx / 2.;
  Scal a;
  if (u <= w / ny) {
    a = std::sqrt(2. * u * nx * ny);
  } else if (u <= 1. - w / ny) {
    a = u * ny + w;
  } else {
    a = nx + ny - std::sqrt(2. * nx * ny * (1. - u));
  }

  if (n[0] < 0.) {
    a += n[0];
  }
  if (n[1] < 0.) {
    a += n[1];
  }

  // shift such that a=0 is cell center
  a -= Vect(0.5).dot(n);

  return a;
}

template <class M>
class Vof : public AdvectionSolver<M> {
  using Mesh = M;
  using P = AdvectionSolver<Mesh>;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  static constexpr size_t dim = Mesh::dim;

  using P::m;
  using P::ffv_;
  using P::fcs_;
  LayersData<FieldCell<Scal>> fc_u_;
  MapFace<std::shared_ptr<ConditionFace>> mf_u_cond_;

  FieldCell<Scal> fc_a_; // alpha (plane constant)
  FieldCell<Vect> fc_n_; // n (normal to plane)
  FieldCell<Scal> fc_us_; // smooth field

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
      , fc_us_(m, 0)
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
            dt * (*fcs_)[c]; // source
      }

      auto p = [&](Vect n, Scal u) {
        std::cout << "n=" << n << " u=" << u << " a=" << GetLineA(n, u)
            << std::endl;
      };
      p(Vect(1., 0., 0.), 0.);
      p(Vect(-1., 0., 0.), 0.);
      p(Vect(1., 0., 0.), 0.5);
      p(Vect(-1., 0., 0.), 0.5);
    }
    const bool split = false;
    for (size_t d = 0; d < (split ? dim : 1); ++d) {
      auto& uc = fc_u_.iter_curr;
      if (sem("pre")) {
        fc_us_ = uc;
      }
      if (sem.Nested("smooth")) {
        Smoothen(fc_us_, mf_u_cond_, m, 1);
      }
      if (sem("adv")) {
        auto us = fc_us_;
        auto uf = Interpolate(us, mf_u_cond_, m);
        auto gc = Gradient(uf, m);

        // commit smooth
        uc = us;

        for (auto c : m.Cells()) {
          const Vect g = gc[c];
          const Vect n =g / (g.norm() + 1e-10);
          fc_n_[c] = n;
          fc_a_[c] = GetLineA(n, uc[c]);
        }

        auto& ffv = *ffv_; // [f]ield [f]ace [v]olume flux
        // ... advection
        m.Comm(&uc);
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
