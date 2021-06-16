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
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Vect3 = generic::Vect<Scal, 3>;
  TaylorGreen() : ModuleInitVelocity<M>("taylor-green") {}
  void operator()(FieldCell<Vect>& fcv, const Vars& var, const M& m) override {
    for (auto i : m.AllCells()) {
      auto& v = fcv[i];
      Vect3 x(m.GetCenter(i));
      if (var.Int["dim"] == 2 && M::dim > 2) {
        x[2] = 0;
      }
      v = Vect(0);
      v[0] = std::sin(x[0]) * std::cos(x[1]) * std::cos(x[2]);
      v[1] = -std::cos(x[0]) * std::sin(x[1]) * std::cos(x[2]);
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
      v = Vect(0);
      v[0] = (x - center).abs() < radius ? 0.5 : -0.5;
      v[1] = std::sin(x[0] * 2 * M_PI / wavelength) * amplitude;
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
      fcv[c] = Vect(0);
      fcv[c][0] = vy0 + (vy1 - vy0) * x[1];
    }
  }
};

template <class M>
class SingleVortex : public ModuleInitVelocity<M> {
 public:
  using Vect = typename M::Vect;
  SingleVortex() : ModuleInitVelocity<M>("single_vortex") {}
  void operator()(FieldCell<Vect>& fcv, const Vars&, const M& m) override {
    using Scal = typename M::Scal;
    for (auto c : m.AllCellsM()) {
      const Scal x = c.center[0] * M_PI;
      const Scal y = c.center[1] * M_PI;
      Scal& u = fcv[c][0];
      Scal& v = fcv[c][1];
      u = std::sin(x) * std::cos(y);
      v = -std::cos(x) * std::sin(y);
    }
  }
};

#define XX(M)                                                             \
  RegisterModule<TaylorGreen<M>>(), RegisterModule<KelvinHelmholtz<M>>(), \
      RegisterModule<Uniform<M>>(), RegisterModule<Couette<M>>(),         \
      RegisterModule<SingleVortex<M>>(),

#define COMMA ,
#define X(dim) XX(MeshCartesian<double COMMA dim>)

bool kReg[] = {MULTIDIMX};

} // namespace init_velocity
