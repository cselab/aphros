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

// Volume fraction field.
// par: parameters
// M: mesh
// Returns:
// std::function<void(GField<Cell>& fc,const M& m)>
// fc: field to fill [i]
// m: mesh
template <class M>
std::function<void(FieldCell<typename M::Scal>&,const M&)> 
CreateInitO(Vars& par, bool verb=true) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  std::function<void(FieldCell<Scal>&,const M&)> g; // result

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

    g = [dim,pp](FieldCell<Scal>& fc, const M& m) { 
      if (pp.empty()) {
        fc.Reinit(m, 0.);
      } else {
        for (auto c : m.Cells()) {
          auto x = m.GetCenter(c);
          for (size_t i = 0; i < pp.size(); ++i) {
            auto& p = pp[i];
            fc[c] = (fm >= 0. ? 1. : 0.);
          }
          // cell size
          Vect h = m.GetNode(m.GetNeighbourNode(c, 7)) - 
              m.GetNode(m.GetNeighbourNode(c, 0));
          // volume fraction
          auto& p = pp[im];
          fc[c] = (fm >= 0. ? 1. : 0.);
        }
      }
    };
  } else if (v == "box") {
    Vect xc(par.Vect["box_c"]);
    Scal s = par.Double["box_s"];
    g = [xc,s](FieldCell<Scal>& fc, const M& m) { 
      for (auto c : m.Cells()) {
        fc[c] = (xc - m.GetCenter(c)).norminf() < s * 0.5 ? 1. : 0.; 
      }
    };
  } else if (v == "zero") {
    g = [](FieldCell<Scal>& fc, const M& m) { 
      for (auto c : m.Cells()) {
        fc[c] = 0.;
      }
    };
  } else {
    if (verb) {
      std::cout << "Notify: no color function for init_vf=" << v << std::endl;
    }
  }
  return g;
}
