#pragma once

#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

#include "func/primlist.h"
#include "geom/block.h"
#include "geom/field.h"
#include "geom/vect.h"
#include "parse/vars.h"

// Init of color function.
// par: parameters
// M: mesh
// Returns:
// std::function<void(GField<Cell>& fc,const M& m)>
// fc: field to fill [i], resized if needed
// m: mesh
template <class M>
std::function<void(
    FieldCell<typename M::Scal>&, const FieldCell<typename M::Scal>&, const M&)>
CreateInitCl(const Vars& par, bool verb = true) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr Scal kClNone = -1.; // no color

  std::string v = par.String["init_vf"];
  if (v == "list") {
    const std::string fn = par.String["list_path"];
    const size_t edim = par.Int["dim"];

    // TODO revise with bcast
    std::ifstream fin(fn);
    if (verb) {
      std::cout << "Open list of primitives '" << fn << "' for colors"
                << std::endl;
    }
    if (!fin.good()) {
      throw std::runtime_error("Can't open list of primitives");
    }
    auto pp = UPrimList<Scal>::Parse(fin, verb, edim);

    return
        [edim, pp](FieldCell<Scal>& cl, const FieldCell<Scal>& vf, const M& m) {
          cl.Reinit(m, kClNone);
          if (pp.empty()) {
            return;
          }

          for (auto c : m.Cells()) {
            if (vf[c] > 0.) {
              auto x = m.GetCenter(c);
              Scal fm = -std::numeric_limits<Scal>::max(); // maximum level-set
              size_t im = 0;
              ; // index of maximum
              for (size_t i = 0; i < pp.size(); ++i) {
                auto& p = pp[i];
                Scal fi = p.ls(x);
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
    return [xc, s](FieldCell<Scal>& fc, const FieldCell<Scal>&, const M& m) {
      fc.Reinit(m, kClNone);
      for (auto c : m.Cells()) {
        if ((xc - m.GetCenter(c)).norminf() < s * 0.5) {
          fc[c] = 1.;
        }
      }
    };
  } else if (v == "grid") {
    using MIdx = typename M::MIdx;
    const int dx = par.Int["initvf_grid_dx"];
    const int dy = par.Int["initvf_grid_dy"];
    const int dz = par.Int["initvf_grid_dz"];
    const MIdx dw(dx, dy, dz);
    return [dw](FieldCell<Scal>& fc, const FieldCell<Scal>&, const M& m) {
      fc.Reinit(m, kClNone);
      const auto& bc = m.GetIndexCells();
      const MIdx gs = m.GetGlobalSize();
      const MIdx ds = (gs + dw - MIdx(1)) / dw;
      for (auto c : m.Cells()) {
        const MIdx w = bc.GetMIdx(c);
        const MIdx rw = w / dw;
        fc[c] = rw[0] + rw[1] * ds[0] + rw[2] * ds[0] * ds[1];
      }
    };
  } else if (v == "zero") {
    return [](FieldCell<Scal>& fc, const FieldCell<Scal>&, const M& m) {
      fc.Reinit(m, kClNone);
    };
  }
  if (verb) {
    std::cout << "Info: no color function for init_vf=" << v << std::endl;
  }
  return std::function<void(
      FieldCell<typename M::Scal>&, const FieldCell<typename M::Scal>&,
      const M&)>();
}
