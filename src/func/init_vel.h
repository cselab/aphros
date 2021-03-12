// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <cmath>
#include <functional>
#include <stdexcept>

#include "geom/field.h"
#include "parse/vars.h"
#include "util/module.h"

// Velocity field.
// par: parameters
// Vect: vector type with operator[] and value_type
// Returns:
// std::function<Scal(Vect x)>
// for pointwise evaluation
template <class Vect, class Scal = typename Vect::value_type>
std::function<Vect(Vect, Scal)> CreateInitVel(const Vars& par) {
  std::function<Vect(Vect, Scal)> f; // result

  std::string v = par.String["init_vel"];
  if (v == "uni") {
    Vect vel(par.Vect["vel"]);
    f = [vel](Vect, Scal) -> Vect { return vel; };
  } else if (v == "sincos") {
    Scal revt = par.Double["revt"]; // reverse time
    f = [revt](Vect x, Scal t) -> Vect {
      x = x * M_PI;
      Vect res(0);
      res[0] = std::sin(x[0]) * std::cos(x[1]);
      res[1] = -std::cos(x[0]) * std::sin(x[1]);
      if (t > revt) {
        res *= -1.;
      }
      return res;
    };
  } else if (v == "stretch") {
    Scal mg = par.Double["stretch_magn"];
    Vect o(par.Vect["stretch_origin"]);
    f = [mg, o](Vect x, Scal) -> Vect {
      x -= o;
      Vect res(0);
      res[0] = x[0];
      res[1] = -x[1];
      return res * mg;
    };
  } else {
    fassert(false, "Unknown init_vel=" + v);
  }

  return f;
}

template <class M>
class ModuleInitVelocity : public Module<ModuleInitVelocity<M>> {
 public:
  using Vect = typename M::Vect;
  using Module<ModuleInitVelocity>::Module;
  virtual void operator()(
      FieldCell<Vect>& fcv, const Vars& var, const M& m) = 0;
};
