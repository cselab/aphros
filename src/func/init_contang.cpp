// Created by Petr Karnakov on 04.05.2020
// Copyright 2020 ETH Zurich

#include <cmath>

#include "geom/mesh.h"
#include "init_contang.h"

// See note about namespaces in module.h

namespace init_contang {

template <class M>
class Uniform : public ModuleInitContang<M> {
 public:
  using Scal = typename M::Scal;
  Uniform() : ModuleInitContang<M>("uniform") {}
  void operator()(
      FieldCell<Scal>& fc_contang, const Vars& var, const M& m) override {
    fc_contang.Reinit(m, var.Double["contang"]);
  }
};

using M = MeshStructured<double, 3>;

bool kReg[] = {
    RegisterModule<Uniform<M>>(),
};

} // namespace init_contang
