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
  using Vect = typename M::Vect;
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
      const Vect x0(conf.x0);
      const Vect x1(conf.x1);
      Scal a = (m.GetCenter(c) - x0).dot(x1 - x0) / (x1 - x0).sqrnorm();
      a = std::max(0., std::min(1., a));
      fc_contang[c] = conf.contang0 * (1 - a) + conf.contang1 * a;
    }
  }
};

template <class M>
class Radial : public ModuleInitContang<M> {
 public:
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
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
      Scal a = m.GetCenter(c).dist(Vect(conf.x0)) / conf.r;
      a = std::min(1., a);
      fc_contang[c] = conf.contang0 * (1 - a) + conf.contang1 * a;
    }
  }
};

#define XX(M)                                                \
  RegisterModule<Uniform<M>>(), RegisterModule<Linear<M>>(), \
      RegisterModule<Radial<M>>(),

#define COMMA ,
#define X(dim) XX(MeshCartesian<double COMMA dim>)

bool kReg[] = {MULTIDIMX};

} // namespace init_contang
