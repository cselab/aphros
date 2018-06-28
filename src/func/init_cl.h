#pragma once

#include <stdexcept>
#include <functional>
#include <cmath>
#include <limits>

#include "parse/vars.h"
#include "geom/field.h"
#include "solver/vof.h"
#include "geom/block.h"
#include "geom/vect.h"
#include "solver/tracker.h"

// Init of color function.
// par: parameters
// M: mesh
// Returns:
// std::function<void(GField<Cell>& fc,const M& m)>
// fc: field to fill [i]
// m: mesh
template <class M>
std::function<void(FieldCell<typename M::Scal>&,
                   const FieldCell<typename M::Scal>&,const M&)> 
CreateInitCl(Vars& par, bool verb=true) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  // result
  std::function<void(FieldCell<Scal>&,const FieldCell<Scal>&,const M&)> g;

  std::string v = par.String["init_vf"];
  if (v == "list") {
    std::string fn = par.String["list_path"];
    size_t dim = par.Int["dim"];

    // TODO: move, same from init_u.h
    struct P { Vect c; Scal r; };
    std::vector<P> pp;

    std::ifstream f(fn);
    if (!f.good()) {
      throw std::runtime_error("Can't open particle list '" + fn + "'");
    }

    // Read until eof
    while (true) {
      P p;
      // Read single particle: x y z r
      f >> p.c[0] >> p.c[1] >> p.c[2] >> p.r;
      if (f.good()) {
        pp.push_back(p);
      } else {
        break;
      }
    }

    if (verb) {
      std::cout << "color: Read " << pp.size() << " particles from " 
          << "'" << fn << "'" << std::endl;
    }

    g = [dim,pp](FieldCell<Scal>& cl, const FieldCell<Scal>& vf, const M& m) { 
      if (pp.empty()) {
        cl.Reinit(m, 0.);
      } else {
        for (auto c : m.Cells()) {
          auto x = m.GetCenter(c);
          // find nearest particle
          // TODO: check if distance to particle edge needed instead
          size_t im = 0;
          for (size_t i = 1; i < pp.size(); ++i) {
            if ((x - pp[i].c).sqrnorm() < (x - pp[im].c).sqrnorm()) {
              im = i;
            }
          }
          // index of nearest particle
          // (assume exact 0 outside particles)
          cl[c] = (vf[c] > 0. ? im : 0.);
        }
      }
    };
  } else if (v == "box") {
    Vect xc(par.Vect["box_c"]);
    Scal s = par.Double["box_s"];
    g = [xc,s](FieldCell<Scal>& fc, const FieldCell<Scal>& fcvf, const M& m) { 
      for (auto c : m.Cells()) {
        fc[c] = (xc - m.GetCenter(c)).norminf() < s * 0.5 ? 1. : 0.; 
      }
    };
  } else if (v == "zero") {
    g = [](FieldCell<Scal>& fc, const FieldCell<Scal>& fcvf, const M& m) { 
      for (auto c : m.Cells()) {
        fc[c] = 0.;
      }
    };
  } else {
    if (verb) {
      std::cout << "Info: no color function for init_vf=" << v << std::endl;
    }
  }
  return g;
}
