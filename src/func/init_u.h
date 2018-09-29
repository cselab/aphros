#pragma once

#include <stdexcept>
#include <functional>
#include <cmath>
#include <limits>
#include <sstream>

#include "parse/vars.h"
#include "geom/field.h"
#include "solver/reconst.h"
#include "geom/block.h"
#include "geom/vect.h"
#include "overlap/over.h"

// Volume fraction cut by interface defined by level-set function
// f: level-set function, interface f=0, f>0 for volume fraction 1
// xc: cell center
// h: cell size
template <class Scal, size_t dim=3>
Scal GetLevelSetVolume(std::function<Scal(const GVect<Scal, 3>&)> f,
                       const GVect<Scal, 3>& xc, const GVect<Scal, 3>& h) {
  using Vect = GVect<Scal, dim>;
  using MIdx = GVect<IntIdx, dim>;

  const Scal dx = 1e-3; // step to compute gradient relative to h

  GBlockCells<dim> b(MIdx(2));

  Scal fc = f(xc);
  for (auto w : b) {
    Vect x = xc + (Vect(w) - Vect(0.5)) * h;
    if ((fc > 0.) != (f(x) > 0.)) {
      // linear approximation
      // f(x) = f(xc) + (x - xc).dot(grad(f, xc))
      Vect n; // normal
      for (size_t i = 0; i < dim; ++i) {
        Vect xp(xc), xm(xc);
        Scal dxh = dx * h[i];
        xp[i] += dxh * 0.5;
        xm[i] -= dxh * 0.5;
        n[i] = (f(xp) - f(xm)) / dxh;
      }
      Scal a; // line constant
      a = f(xc);

      return Reconst<Scal>::GetLineU(n, a, h);
    }
  }

  return fc > 0. ? 1. : 0.;
}

// Volume fraction field.
// par: parameters
// M: mesh
// Returns:
// std::function<void(GField<Cell>& fc,const M& m)>
// fc: field to fill [i]
// m: mesh
template <class M>
std::function<void(FieldCell<typename M::Scal>&,const M&)> 
CreateInitU(Vars& par, bool verb=true) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  std::function<void(FieldCell<Scal>&,const M&)> g; // result

  std::string v = par.String["init_vf"];
  if (v == "circle") {
    Vect xc = Vect(par.Vect["circle_c"]);
    Scal r = par.Double["circle_r"];
    size_t dim = par.Int["dim"];
    g = [xc,r,dim](FieldCell<Scal>& fc, const M& m) { 
      for (auto c : m.Cells()) {
        auto dx = m.GetCenter(c) - xc;
        if (dim == 2) {
          dx[2] = 0.;
        }
        fc[c] = dx.sqrnorm() < sqr(r) ? 1. : 0.;
      }
    };
  } else if (v == "circlels") {
    Vect xc = Vect(par.Vect["circle_c"]);
    Scal r = par.Double["circle_r"];
    size_t dim = par.Int["dim"];
    g = [xc,r,dim](FieldCell<Scal>& fc, const M& m) { 
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
      for (auto c : m.Cells()) {
        auto x = m.GetCenter(c);
        Vect h = m.GetNode(m.GetNeighbourNode(c, 7)) - 
            m.GetNode(m.GetNeighbourNode(c, 0));
        fc[c] = GetLevelSetVolume<Scal>(f, (x - xc) / r, h / r);
      }
    };
  } else if (v == "list") {
    std::string fn = par.String["list_path"];
    int ls = par.Int["list_ls"]; // 0: stepwise, 1: level-set, 2: overlap
    size_t dim = par.Int["dim"];

    // elliptic partilces TODO: generalize
    struct P { 
      Vect c;  // center
      Vect r;  // axes in coordinate directions
    };
    std::vector<P> pp;

    std::ifstream f(fn);
    if (!f.good()) {
      throw std::runtime_error("Can't open particle list '" + fn + "'");
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
        s = c + s; // append first character
        std::stringstream st(s);
        st >> p.c[0] >> p.c[1] >> p.c[2];
        st >> p.r[0];
        if (st.fail()) {
          throw std::runtime_error("list: missing rx in '" + s + "'");
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
      std::cout << "Read " << pp.size() << " particles from " 
          << "'" << fn << "'" << std::endl;
    }

    g = [dim,ls,pp](FieldCell<Scal>& fc, const M& m) { 
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
      if (pp.empty()) {
        fc.Reinit(m, 0.);
      } else {
        for (auto c : m.Cells()) {
          auto x = m.GetCenter(c);
          // find particle with largest value of level-set
          Scal fm = -std::numeric_limits<Scal>::max(); 
          size_t im; // index of particle
          for (size_t i = 0; i < pp.size(); ++i) {
            auto& p = pp[i];
            Scal fi = f((x - p.c) / p.r);
            if (fi > fm) {
              fm = fi;
              im = i;
            }
          }
          // cell size
          Vect h = m.GetNode(m.GetNeighbourNode(c, 7)) - 
              m.GetNode(m.GetNeighbourNode(c, 0));
          // volume fraction
          auto& p = pp[im];
          if (ls == 1) {
            fc[c] = GetLevelSetVolume<Scal>(f, (x - p.c) / p.r, h / p.r);
          } else if (ls == 2) {
            fc[c] = GetSphereOverlap((x - p.c) / p.r, h / p.r, Vect(0), 1.);
          } else {
            fc[c] = (fm >= 0. ? 1. : 0.);
          }
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
  } else if (v == "line") {
    Vect xc(par.Vect["line_c"]); // center
    Vect n(par.Vect["line_n"]);  // normal
    Scal h(par.Double["line_h"]);  // thickness 
    g = [xc,n,h](FieldCell<Scal>& fc, const M& m) { 
      for (auto c : m.Cells()) {
        Scal d = (m.GetCenter(c) - xc).dot(n);
        fc[c] = 1. / (1. + std::exp(-d / h));
      }
    };
  } else if (v == "sin") {
    Vect k;
    if (auto p = par.Vect("sin_k")) {
      k = Vect(*p);
    } else {
      k = Vect(2. * M_PI);
    }

    g = [k](FieldCell<Scal>& fc, const M& m) { 
      for (auto c : m.Cells()) {
        Vect z = m.GetCenter(c) * k;
        fc[c] = std::sin(z[0]) * std::sin(z[1]) * std::sin(z[2]);
      }
    };
  } else if (v == "sinc") {
    Vect k(par.Vect["sinc_k"]);
    g = [k](FieldCell<Scal>& fc, const M& m) { 
      for (auto c : m.Cells()) {
        Vect x = m.GetCenter(c);
        x -= Vect(0.5);
        x *= k;
        Scal r = x.norm();
        Scal u0 = -0.2;
        Scal u1 = 1.;
        Scal u = std::sin(r) / r;
        u = (u - u0) / (u1 - u0);
        fc[c] = std::max(0., std::min(1., u));
      }
    };
  } else if (v == "zero") {
    g = [](FieldCell<Scal>& fc, const M& m) { 
      for (auto c : m.Cells()) {
        fc[c] = 0.;
      }
    };
  } else {
    throw std::runtime_error("Unknown init_vf=" + v);
  }
  return g;
}
