#pragma once

#include <stdexcept>
#include <functional>
#include <cmath>

#include "parse/vars.h"
#include "geom/field.h"

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
    Vect c;
    Scal r;
    if (par.Vect("circle_c") && par.Double("circle_r")) {
      c = Vect(par.Vect["circle_c"]);
      r = par.Double["circle_r"];
    } else {
      c = Vect(0.5, 0.263662, 0.);
      r = 0.2;
    }
    f = [c,r](Vect x) -> Scal { 
      return c.dist(x) < r ? 1. : 0.; 
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

// Volume fraction field.
// par: parameters
// M: mesh
// Returns:
// std::function<void(GField<Cell>& fc)>
template <class M>
std::function<void(FieldCell<typename M::Scal>&, const M&)> 
CreateInitU(Vars& par) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  std::function<void(FieldCell<Scal>&,M&)> g; // result

  std::string v = par.String["init_vf"];
  if (v == "circle") {
    Vect xc = Vect(par.Vect["circle_c"]);
    Scal r = par.Double["circle_r"];
    g = [xc,r](FieldCell<Scal>& fc, const M& m) -> Scal { 
      for (auto c : m.Cells()) {
        fc[c] = (xc.dist(m.GetCenter(c)) < r ? 1. : 0.);
      }
    };
  } else {
    throw std::runtime_error("Unknown init_vf=" + v);
  }
  return g;
}
