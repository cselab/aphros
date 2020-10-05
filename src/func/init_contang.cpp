// Created by Petr Karnakov on 04.05.2020
// Copyright 2020 ETH Zurich

#include <cmath>

#include "geom/mesh.h"
#include "init_contang.h"
#include "parse/config.h"

DECLARE_FORCE_LINK_TARGET(init_contang);

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
    for (auto c : m.AllCells()) {
      Scal a = (m.GetCenter(c) - conf.x0).dot(conf.x1 - conf.x0) /
                     (conf.x1 - conf.x0).sqrnorm();
      a = std::max(0., std::min(1., a));
      fc_contang[c] = conf.contang0 * (1 - a) + conf.contang1 * a;
    }
  }
};

template <class M>
class Radial : public ModuleInitContang<M> {
 public:
  using Scal = typename M::Scal;
  Radial() : ModuleInitContang<M>("radial") {}
  void operator()(
      FieldCell<Scal>& fc_contang, const Vars& var, const M& m) override {
    struct : ConfigBase {
      VAR_VECT3(x0);
      VAR_DOUBLE(r);
      VAR_DOUBLE(contang0);
      VAR_DOUBLE(contang1);
    } conf;
    conf.Read(var, "contang_");
    conf.contang0 *= M_PI / 180.;
    conf.contang1 *= M_PI / 180.;
    for (auto c : m.AllCells()) {
      Scal a = m.GetCenter(c).dist(conf.x0) / conf.r;
      a = std::min(1., a);
      fc_contang[c] = conf.contang0 * (1 - a) + conf.contang1 * a;
    }
  }
};

using M = MeshStructured<double, 3>;

bool kReg[] = {
    RegisterModule<Uniform<M>>(),
    RegisterModule<Linear<M>>(),
    RegisterModule<Radial<M>>(),
};

} // namespace init_contang
