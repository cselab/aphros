#pragma once

#include <stdexcept>
#include <functional>
#include <cmath>

#include "hydro/field.h"
#include "CubismDistr/Vars.h"



// Volume fraction field.
// par: parameters
// Vect: vector type with operator[] and value_type
// Returns:
// std::function<Scal(Vect x)>
// for pointwise evaluation
template <class Vect, class Scal=typename Vect::value_type>
std::function<Scal(Vect)> CreateInitU(Vars& par) {
  std::function<Scal(Vect)> f; // result

  std::string v = par.String["init_vf"];
  if (v == "circle") {
    f = [](Vect x) -> Scal { 
      return Vect(0.5, 0.263662, 0.).dist(x) < 0.2 ? 1. : 0.; 
    };
  } else if (v == "box") {
    Vect c(par.Vect["box_center"]);
    Scal s = par.Double["box_size"];
    f = [c,s](Vect x) -> Scal { 
      return (x - c).norminf() < s * 0.5 ? 1. : 0.; 
    };
  } else if (v == "line") {
    Vect xc(par.Vect["linec"]); // center
    Vect n(par.Vect["linen"]);  // normal
    Scal s(par.Double["lines"]);  // sigma (sharpness)
    f = [xc, n, s](Vect x) -> Scal { 
      Scal u = (x - xc).dot(n);
      u = 1. / (1. + std::exp(-s * u));
      return u;
    };
  } else if (v == "sinc") {
    Vect k(par.Vect["sinck"]);
    f = [k](Vect x) -> Scal { 
      x -= Vect(0.5);
      x *= k;
      Scal r = x.norm();
      Scal u0 = -0.2;
      Scal u1 = 1.;
      Scal u = std::sin(r) / r;
      u = (u - u0) / (u1 - u0);
      return std::max(0., std::min(1., u));
    };
  } else {
    throw std::runtime_error("Unknown init_vf=" + v);
  }

  return f;
}
