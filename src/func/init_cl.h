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
// fc: field to fill [i], resized if needed
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

    // elliptic partilces TODO: generalize
    struct P { 
      Vect c;  // center
      Vect r;  // axes in coordinate directions
    };
    std::vector<P> pp;

    std::ifstream f(fn);
    if (!f.good()) {
      throw std::runtime_error("color: Can't open particle list '" + fn + "'");
    }

    f >> std::skipws;
    // Read until eof
    while (true) {
      P p;
      // Read single particle: x y z r
      // first character (to skip empty strings)
      char c;
      f >> c;
      if (f.good()) {
        // still have lines
        std::string s;
        std::getline(f, s);
        if (c == '#') {
          continue;
        }
        s = c + s; // append first character
        std::stringstream st(s);
        st >> p.c[0] >> p.c[1] >> p.c[2];
        st >> p.r[0];
        if (st.fail()) {
          throw std::runtime_error("color_list: missing rx in '" + s + "'");
        }
        st >> p.r[1];
        if (st.fail()) {
          p.r[1] = p.r[0];
          p.r[2] = p.r[0];
        } else {
          st >> p.r[2];
          if (st.fail()) {
            p.r[2] = p.r[0];
          }
        }
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
      // level-set for particle of radius 1 centered at zero,
      // positive inside,
      // cylinder along z if dim=2
      auto f = [dim](const Vect& x) -> Scal {
        Vect xd = x;
        if (dim == 2) {
          xd[2] = 0.;
        }
        return 1. - xd.sqrnorm();
      };

      Scal kNone = solver::Tracker<M>::kNone;
      cl.Reinit(m, kNone);
      if (pp.empty()) {
        return;
      }

      for (auto c : m.Cells()) {
        if (vf[c] > 0.) { // liquid in cell // TODO: threshold
          auto x = m.GetCenter(c);
          // find particle with largest value of level-set
          Scal fm = -std::numeric_limits<Scal>::max(); 
          size_t im = 0;; // index of particle
          for (size_t i = 0; i < pp.size(); ++i) {
            auto& p = pp[i];
            Scal fi = f((x - p.c) / p.r);
            if (fi > fm) {
              fm = fi;
              im = i;
            }
          }
          cl[c] = im;
        }
      }
    };
  } else if (v == "box") {
    Vect xc(par.Vect["box_c"]);
    Scal s = par.Double["box_s"];
    g = [xc,s](FieldCell<Scal>& fc, const FieldCell<Scal>& /*fcvf*/, const M& m) { 
      Scal kNone = solver::Tracker<M>::kNone;
      fc.Reinit(m, kNone);
      for (auto c : m.Cells()) {
        if ((xc - m.GetCenter(c)).norminf() < s * 0.5) {
          fc[c] = 1.; 
        }
      }
    };
  } else if (v == "zero") {
    g = [](FieldCell<Scal>& fc, const FieldCell<Scal>& /*fcvf*/, const M& m) { 
      Scal kNone = solver::Tracker<M>::kNone;
      fc.Reinit(m, kNone);
    };
  } else {
    if (verb) {
      std::cout << "Info: no color function for init_vf=" << v << std::endl;
    }
  }
  return g;
}
