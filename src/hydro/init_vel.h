#pragma once

#include <stdexcept>
#include <functional>
#include <cmath>

#include "CubismDistr/Vars.h"

template <class M>
using FuncVVS = 
    std::function<typename M::Vect(typename M::Vect, typename M::Scal)>;

template <class M>
FuncVVS<M> CreateInitVel(Vars& par) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  FuncVVS<M> f; // result

  std::string v = par.String["init_vel"];
  if (v == "uni") {
    Vect vel(par.Vect["vel"]);
    f = [vel](Vect x, Scal /*t*/) -> Vect { 
      return Vect(vel); 
    };
  } else if (v == "sincos") {
    Scal revt = par.Double["revt"]; // reverse time
    f = [revt](Vect x, Scal t) -> Vect { 
      x = x * M_PI;
      Vect r(std::sin(x[0]) * std::cos(x[1]), 
             -std::cos(x[0]) * std::sin(x[1]),
             0.);
      if (t > revt) { r *= -1.; }
      return r; 
    };
  } else if (v == "stretch") {
    Scal mg = par.Double["stretch_magn"];
    Vect o(par.Vect["stretch_origin"]);
    f = [mg, o](Vect x, Scal t) -> Vect { 
      x -= o;
      return Vect(x[0], -x[1], 0.) * mg;
    };
  } else {
    throw std::runtime_error("Unknown init_vel=" + v);
  }

  return f;
}
