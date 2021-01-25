// Created by Petr Karnakov on 20.04.2020
// Copyright 2020 ETH Zurich

#include <cmath>

#include "geom/mesh.h"
#include "init_vel.h"

DECLARE_FORCE_LINK_TARGET(init_vel);

// See note about namespaces in module.h
namespace init_velocity {

template <class M>
class TaylorGreen : public ModuleInitVelocity<M> {
 public:
  using Vect = typename M::Vect;
  TaylorGreen() : ModuleInitVelocity<M>("taylor-green") {}
  void operator()(FieldCell<Vect>& fcv, const Vars& var, const M& m) override {
    for (auto i : m.AllCells()) {
      auto& v = fcv[i];
      auto x = m.GetCenter(i);
      if (var.Int["dim"] == 2) {
        x[2] = 0.;
      }
      v[0] = std::sin(x[0]) * std::cos(x[1]) * std::cos(x[2]);
      v[1] = -std::cos(x[0]) * std::sin(x[1]) * std::cos(x[2]);
      v[2] = 0.;
    }
  }
};

template <class M>
class KelvinHelmholtz : public ModuleInitVelocity<M> {
 public:
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  KelvinHelmholtz() : ModuleInitVelocity<M>("kelvin-helmholtz") {}
  void operator()(FieldCell<Vect>& fcv, const Vars& var, const M& m) override {
    const Vect center(var.Vect["kh_center"]);
    const Vect radius(var.Vect["kh_radius"]);
    const Scal wavelength = var.Double["kh_wavelength"];
    const Scal amplitude = var.Double["kh_amplitude"];
    for (auto i : m.AllCells()) {
      Vect& v = fcv[i];
      const Vect x = m.GetCenter(i);
      v[0] = (x - center).abs() < radius ? 0.5 : -0.5;
      v[1] = std::sin(x[0] * 2 * M_PI / wavelength) * amplitude;
      v[2] = 0.;
    }
  }
};

template <class M>
class Uniform : public ModuleInitVelocity<M> {
 public:
  using Vect = typename M::Vect;
  Uniform() : ModuleInitVelocity<M>("uniform") {}
  void operator()(FieldCell<Vect>& fcv, const Vars& var, const M& m) override {
    const Vect v(var.Vect["vel"]);
    for (auto c : m.AllCells()) {
      fcv[c] = v;
    }
  }
};

template <class M>
class Couette : public ModuleInitVelocity<M> {
 public:
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  Couette() : ModuleInitVelocity<M>("couette") {}
  void operator()(FieldCell<Vect>& fcv, const Vars& var, const M& m) override {
    const Scal vy0 = var.Double["vx_y0"];
    const Scal vy1 = var.Double["vx_y1"];
    for (auto c : m.AllCells()) {
      const Vect x = m.GetCenter(c);
      fcv[c] = Vect(vy0 + (vy1 - vy0) * x[1], 0., 0.);
    }
  }
};

using M = MeshStructured<double, 3>;

bool kReg[] = {
    RegisterModule<TaylorGreen<M>>(),
    RegisterModule<KelvinHelmholtz<M>>(),
    RegisterModule<Uniform<M>>(),
    RegisterModule<Couette<M>>(),
};

} // namespace init_velocity
