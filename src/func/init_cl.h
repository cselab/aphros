#pragma once

#include <stdexcept>
#include <functional>
#include <cmath>
#include <limits>

#include "parse/vars.h"
#include "geom/field.h"
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

  std::string v = par.String["init_vf"];
  if (v == "list") {
    std::string fn = par.String["list_path"];
    size_t dim = par.Int["dim"];

    auto pp = UPrimList<Scal>::Parse(par.String["list_path"], verb);

    return [dim,pp](FieldCell<Scal>& cl, const FieldCell<Scal>& vf, 
                    const M& m) { 
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
    return [xc,s](FieldCell<Scal>& fc, const FieldCell<Scal>&, const M& m) {
      Scal kNone = solver::Tracker<M>::kNone;
      fc.Reinit(m, kNone);
      for (auto c : m.Cells()) {
        if ((xc - m.GetCenter(c)).norminf() < s * 0.5) {
          fc[c] = 1.; 
        }
      }
    };
  } else if (v == "zero") {
    return [](FieldCell<Scal>& fc, const FieldCell<Scal>&, const M& m) { 
      Scal kNone = solver::Tracker<M>::kNone;
      fc.Reinit(m, kNone);
    };
  }
  if (verb) {
    std::cout << "Info: no color function for init_vf=" << v << std::endl;
  }
  return std::function<void(FieldCell<typename M::Scal>&, 
                            const FieldCell<typename M::Scal>&, 
                            const M&)>();
}
