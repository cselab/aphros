// Created by Petr Karnakov on 04.05.2020
// Copyright 2020 ETH Zurich

#include <cmath>

#include "geom/mesh.h"
#include "init_contang.h"
#include "parse/config.h"

// See note about namespaces in module.h

namespace init_contang {

template <class M>
class Uniform : public ModuleInitContang<M> {
 public:
  using Scal = typename M::Scal;
  Uniform() : ModuleInitContang<M>("uniform") {}
  void operator()(
      FieldCell<Scal>& fc_contang, const Vars& var, const M& m) override {
    fc_contang.Reinit(m, var.Double["contang"] * M_PI / 180);
  }
};

template <class M>
class Linear : public ModuleInitContang<M> {
 public:
  using Scal = typename M::Scal;
  Linear() : ModuleInitContang<M>("linear") {}
  void operator()(
      FieldCell<Scal>& fc_contang, const Vars& var, const M& m) override {
    struct : ConfigBase {
      VAR_VECT3(x0);
      VAR_VECT3(x1);
      VAR_DOUBLE(contang0);
      VAR_DOUBLE(contang1);
    } conf;
    conf.Read(var, "contang_");
    conf.contang0 *= M_PI / 180.;
    conf.contang1 *= M_PI / 180.;
    using Scal = typename M::Scal;
    for (auto c : m.AllCells()) {
      const Scal a = (conf.x1 - m.GetCenter(c)).dot(conf.x1 - conf.x0) /
                     (conf.x1 - conf.x0).sqrnorm();
      fc_contang[c] = conf.contang0 * a + conf.contang1 * (1 - a);
    }
  }
};

using M = MeshStructured<double, 3>;

bool kReg[] = {
    RegisterModule<Uniform<M>>(),
    RegisterModule<Linear<M>>(),
};

} // namespace init_contang
