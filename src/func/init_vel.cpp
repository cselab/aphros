// Created by Petr Karnakov on 20.04.2020
// Copyright 2020 ETH Zurich

#include <cmath>

#include "init_vel.h"
#include "geom/mesh.h"

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
class Uniform : public ModuleInitVelocity<M> {
 public:
  using Vect = typename M::Vect;
  Uniform() : ModuleInitVelocity<M>("uniform") {}
  void operator()(FieldCell<Vect>& fcv, const Vars& var, const M& m) override {
    Vect v(var.Vect["vel"]);
    for (auto c : m.AllCells()) {
      fcv[c] = v;
    }
  }
};

using M = MeshStructured<double, 3>;

bool kReg[] = {
    RegisterModule<TaylorGreen<M>>(),
    RegisterModule<Uniform<M>>(),
};

} // namespace init_velocity
