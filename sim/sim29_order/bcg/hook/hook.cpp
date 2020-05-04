// Created by Petr Karnakov on 04.03.2020
// Copyright 2020 ETH Zurich

#include <iostream>

#include <util/posthook.h>
#include <util/module.h>
#include <func/init_vel.h>

template <class M>
class Bcg : public ModuleInitVelocity<M> {
 public:
  using Vect = typename M::Vect;
  Bcg() : ModuleInitVelocity<M>("bcg") {}
  void operator()(FieldCell<Vect>& fcv, const Vars&, const M& m) override {
    using Scal = typename M::Scal;
    for (auto c : m.AllCells()) {
      // psi = 1/pi * sin**2(pi*x) * sin**2(pi*y)
      const Vect xx = m.GetCenter(c);
      const Scal x = xx[0];
      const Scal y = xx[1];
      Scal& u = fcv[c][0];
      Scal& v = fcv[c][1];
      const Scal pi = M_PI;
      using std::sin;
      using std::cos;
      u = -2 * sqr(sin(pi * x)) * sin(pi * y) * cos(pi * y);
      v = 2 * sqr(sin(pi * y)) * sin(pi * x) * cos(pi * x);
    }
  }
};

using M = MeshStructured<double, 3>;

bool reg[] = {
    RegisterModule<Bcg<M>>(),
};
