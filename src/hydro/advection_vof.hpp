#pragma once

#include <exception>
#include <fstream>

#include "mesh.hpp"
#include "linear.hpp"
#include "solver.hpp"
#include "advection.hpp"

namespace solver {

using namespace geom;

template <class Scal>
inline void Clip(Scal& a, Scal l, Scal u) {
  a = std::max(l, std::min(u, a));
}

template <class Scal>
inline void Clip(Scal& a) {
  Clip(a, 0., 1.);
};

// GetLineA() helper
// assuming 0 < u < 0.5, 0 < nx < ny
template <class Scal>
inline Scal GetLineA0(Scal nx, Scal ny, Scal u) {
  Scal u1 = 0.5 * nx / ny;
  if (u <= u1) {
    return -0.5 * (nx + ny) + std::sqrt(2. * nx * ny * u);
  } else {
    return ny * (u - 0.5);
  }
}

// Volume fraction to line constant.
// n : unit normal
// u: volume fraction
// Returns:
// a: line constant
// Equation of reconstructed plane would be
// ((x-xc)/h).dot(n) = a
// where 
// xc: cell center
// h: cell size
template <class Scal>
inline Scal GetLineA(const GVect<Scal, 3>& n, Scal u) {
  using Vect = GVect<Scal, 3>;

  Scal nx = std::abs(n[0]);
  Scal ny = std::abs(n[1]);
  if (ny < nx) {
    std::swap(nx, ny);
  }
  
  Clip(u);

  if (u < 0.5) {
    return GetLineA0(nx, ny, u);
  } else {
    return -GetLineA0(nx, ny, 1. - u);
  }
}

// GetLineU() helper
// assuming a < 0, 0 < nx < ny
template <class Scal>
inline Scal GetLineU0(Scal nx, Scal ny, Scal a) {
  Scal a1 = 0.5 * (nx - ny);
  if (a <= a1) {
    if (nx == 0.) { // TODO revise
      return 0.;
    } else {
      return std::pow(a + 0.5 * (nx + ny), 2) / (2. * nx * ny);
    }
  } else {
    return 0.5 + a / ny;
  }
}

// Line constant to volume fraction
// n : unit normal
// a: line constant
// Returns:
// u: volume fraction
// (see equation in GetLineA())
template <class Scal>
inline Scal GetLineU(const GVect<Scal, 3>& n, Scal a) {
  using Vect = GVect<Scal, 3>;

  Scal nx = std::abs(n[0]);
  Scal ny = std::abs(n[1]);
  if (ny < nx) {
    std::swap(nx, ny);
  }
  
  Clip(a, -0.5 * (nx + ny), 0.5 * (nx + ny));

  if (a < 0.) {
    return GetLineU0(nx, ny, a);
  } else {
    return 1. - GetLineU0(nx, ny, -a);
  }
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
        Scal a = GetLineA(n, u);
        std::cout 
            << "n=" << n << " u=" << u 
            << " a=" << a
            << " u(a)=" << GetLineU(n, a) 
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
