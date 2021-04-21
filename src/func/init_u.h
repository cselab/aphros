// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <cmath>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "debug/isnan.h"
#include "dump/hdf.h"
#include "geom/block.h"
#include "geom/field.h"
#include "geom/vect.h"
#include "parse/vars.h"
#include "primlist.h"
#include "solver/reconst.h"

// Volume fraction cut by interface defined by level-set function
// ls: level-set function, interface ls=0, ls>0 for volume fraction 1
// xc: cell center
// h: cell size
template <class Scal, size_t dim>
static Scal GetLevelSetVolume(
    std::function<Scal(const generic::Vect<Scal, dim>&)> ls,
    const generic::Vect<Scal, dim>& xc, const generic::Vect<Scal, dim>& h) {
  using Vect = generic::Vect<Scal, dim>;
  using MIdx = generic::Vect<IntIdx, dim>;

  const Scal dx = 1e-3; // step to compute gradient relative to h

  GBlock<size_t, dim> b(MIdx(2));
  const Scal lsc = ls(xc);
  // traverse nodes of cell
  for (auto w : b) {
    const Vect x = xc + (Vect(w) - Vect(0.5)) * h;
    if ((lsc > 0.) != (ls(x) > 0.)) { // cell crossed by interface
      // linear approximation
      // ls(x) = ls(xc) + (x - xc).dot(grad(ls, xc))
      Vect n; // normal
      for (size_t i = 0; i < dim; ++i) {
        Vect xp(xc), xm(xc);
        const Scal dxh = dx * h[i];
        xp[i] += dxh * 0.5;
        xm[i] -= dxh * 0.5;
        n[i] = (ls(xp) - ls(xm)) / dxh;
      }
      return Reconst<Scal>::GetLineU(n, lsc, h);
    }
  }

  return lsc > 0. ? 1. : 0.;
}

// Returns point at which interpolant has value 0.
// x0,x1: points
// f0,f1: values
template <class Scal, class Vect>
Vect GetIso(Vect x0, Vect x1, Scal f0, Scal f1) {
  return (x0 * f1 - x1 * f0) / (f1 - f0);
}

template <class Scal>
Scal GetPositiveEdgeFraction(Scal l0, Scal l1) {
  if (l0 * l1 < 0) {
    return l0 < l1 ? l1 / (l1 - l0) : l0 / (l0 - l1);
  }
  return l0 + l1 > 0;
}

// Returns fraction of face area for which ls>0.
// ee: fraction of edge length for which ls>0
//       3
//   -------
//  0|     |
//   |     | 2
//   -------
//     1
template <class Scal>
Scal GetFaceAreaFraction(std::array<Scal, 4> ee) {
  using R = Reconst<Scal>;
  // face center is x=0 y=0
  // line equation
  //   n.dot(x) = a
  if (ee[0] < ee[2]) std::swap(ee[0], ee[2]);
  if (ee[1] < ee[3]) std::swap(ee[1], ee[3]);

  constexpr size_t ni = 4;

  const std::array<Scal, ni> xx{-0.5, ee[1] - 0.5, 0.5, ee[3] - 0.5};
  const std::array<Scal, ni> yy{ee[0] - 0.5, -0.5, ee[2] - 0.5, 0.5};

  // normal
  Scal nx = ee[0] - ee[2];
  Scal ny = ee[1] - ee[3];

  Scal a = 0;
  size_t aw = 0;
  for (size_t i = 0; i < ni; ++i) {
    const auto& e = ee[i];
    if (e > 0 && e < 1) {
      a += xx[i] * nx + yy[i] * ny;
      ++aw;
    }
  }
  if (aw == 0 || nx + ny == 0) {
    return (ee[0] + ee[1] + ee[2] + ee[3]) / ni;
  }
  a /= aw;

  if (nx > ny) {
    std::swap(nx, ny);
  }
  using Vect2 = generic::Vect<Scal, 2>;
  if (a < 0) {
    return R::GetLineU0(Vect2{nx, ny}, a);
  }
  return 1 - R::GetLineU0(Vect2{nx, ny}, -a);
}

// Copmutes face area for which ls > 0.
// fnl: level-set function ls on nodes [i]
template <class M>
FieldFace<typename M::Scal> GetPositiveAreaFraction(
    const FieldNode<typename M::Scal>& fnl, const M& m) {
  using Scal = typename M::Scal;
  FieldFace<Scal> ffs(m, 0);
  for (auto f : m.Faces()) {
    constexpr size_t em = 4;
    std::array<Scal, em> ll;
    for (size_t e = 0; e < em; ++e) {
      const size_t ep = (e + 1) % em;
      const IdxNode n = m.GetNode(f, e);
      const IdxNode np = m.GetNode(f, ep);
      ll[e] = GetPositiveEdgeFraction(fnl[n], fnl[np]);
      ffs[f] = GetFaceAreaFraction(ll);
    }
  }
  return ffs;
}

// fnl: level-set function ls computed on nodes [i]
template <class M>
FieldCell<typename M::Scal> GetPositiveVolumeFraction(
    const FieldNode<typename M::Scal>& fnl, const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using R = Reconst<Scal>;
  auto ffs = GetPositiveAreaFraction(fnl, m);
  FieldCell<Scal> fcu(m, 0);
  for (auto c : m.Cells()) {
    Vect nn(0); // normal
    for (auto q : m.Nci(c)) {
      const IdxFace f = m.GetFace(c, q);
      nn += m.GetNormal(f) * ffs[f] * m.GetOutwardFactor(c, q);
    }
    Scal u = 0;
    if (nn.norm1() == 0) {
      for (auto q : m.Nci(c)) {
        u += ffs[m.GetFace(c, q)];
      }
      u /= m.Nci(c).size();
    } else {
      nn /= -nn.norm1();
      Scal a = 0; // plane constant
      size_t aw = 0;
      for (auto q : m.Nci(c)) {
        const IdxFace f = m.GetFace(c, q);
        const size_t em = m.GetNumNodes(f);
        for (size_t e = 0; e < em; ++e) {
          const size_t ep = (e + 1) % em; // next node
          const IdxNode n = m.GetNode(f, e);
          const IdxNode np = m.GetNode(f, ep);
          const Scal l = fnl[n];
          const Scal lp = fnl[np];
          const Vect x = m.GetNode(n);
          const Vect xp = m.GetNode(np);

          if (l * lp < 0) {
            a += nn.dot(GetIso(x, xp, l, lp) - m.GetCenter(c));
            aw += 1;
          }
        }
      }
      assert(aw > 0);
      a /= aw;
      u = R::GetLineU(nn, a, m.GetCellSize());
    }
    fcu[c] = u;
  }
  return fcu;
}

// Fills volume fraction field from list of primitives.
// fc: field to fill
// list: list of primitives
// edim: effective dimension
// approx: 0: stepwise, 1: level-set, 2: overlap, 3: level-set on nodes
template <class M>
void InitVfList(
    FieldCell<typename M::Scal>& fc, std::istream& list, int approx,
    size_t edim, const M& m, bool verbose);

// Volume fraction field.
// par: parameters
// M: mesh
// Returns:
// std::function<void(GField<Cell>& fc,const M& m)>
// fc: field to fill [i]
// m: mesh
template <class M>
std::function<void(FieldCell<typename M::Scal>&, const M&)> CreateInitU(
    const Vars& par, bool verb) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  (void)verb;

  std::string v = par.String["init_vf"];
  if (v == "circle") {
    const Vect xc(par.Vect["circle_c"]);
    const Scal r = par.Double["circle_r"];
    size_t dim = par.Int["dim"];
    return [xc, r, dim](FieldCell<Scal>& fc, const M& m) {
      for (auto c : m.Cells()) {
        auto dx = m.GetCenter(c) - xc;
        if (dim == 2) {
          dx[2] = 0;
        }
        fc[c] = dx.sqrnorm() < sqr(r) ? 1. : 0.;
      }
    };
  } else if (v == "radial_trapezoid") {
    const Vect xc(par.Vect[v + "_c"]);
    const Scal rmin = par.Double[v + "_rmin"];
    const Scal rmax = par.Double[v + "_rmax"];
    const size_t dim = par.Int["dim"];
    return [xc, rmin, rmax, dim](FieldCell<Scal>& fc, const M& m) {
      for (auto c : m.Cells()) {
        auto dx = m.GetCenter(c) - xc;
        if (M::dim > 2 && dim == 2) {
          dx[2] = 0;
        }
        const Scal r = dx.norm();
        fc[c] = std::max(0., std::min(1., (rmax - r) / (rmax - rmin)));
      }
    };
  } else if (v == "circlels") {
    const Vect xc(par.Vect["circle_c"]);
    const Scal r = par.Double["circle_r"];
    const size_t dim = par.Int["dim"];
    return [xc, r, dim](FieldCell<Scal>& fc, const M& m) {
      // level-set for particle of radius 1 centered at zero,
      // positive inside,
      // cylinder along z if dim=2
      auto f = [dim](const Vect& x) -> Scal {
        auto xd = x;
        if (dim == 2) {
          xd[2] = 0.;
        }
        return 1 - xd.sqrnorm();
      };
      for (auto c : m.Cells()) {
        auto x = m.GetCenter(c);
        Vect h = m.GetCellSize();
        fc[c] = GetLevelSetVolume<Scal, M::dim>(f, (x - xc) / r, h / r);
      }
    };
  } else if (v == "gear") {
    const Vect xc(par.Vect["gear_c"]);
    const Scal r = par.Double["gear_r"];
    const Scal a = par.Double["gear_amp"]; // amplitude relative to r
    const Scal n = par.Double["gear_n"]; // number of peaks
    const size_t dim = par.Int["dim"];
    return [xc, r, dim, a, n](FieldCell<Scal>& fc, const M& m) {
      // level-set for particle of radius 1 centered at zero,
      // modulated by sine-wave of amplitude a and number of periods n
      // positive inside,
      // cylinder along z if dim=2
      auto f = [dim, n, a](const Vect& x) -> Scal {
        Vect xd = x;
        if (M::dim > 2 && dim == 2) {
          xd[2] = 0;
        }
        const Scal phi = std::atan2(xd[1], xd[0]);
        return 1 + std::sin(phi * n) * a - xd.norm();
      };
      for (auto c : m.Cells()) {
        auto x = m.GetCenter(c);
        Vect h = m.GetCellSize();
        fc[c] = GetLevelSetVolume<Scal, M::dim>(f, (x - xc) / r, h / r);
      }
    };
  } else if (v == "soliton") {
    const Scal xc(par.Double["soliton_xc"]);
    const Scal yc(par.Double["soliton_yc"]);
    const Scal yh(par.Double["soliton_yh"]);
    return [xc, yc, yh](FieldCell<Scal>& fc, const M& m) {
      auto f = [xc, yh, yc](const Vect& xx) -> Scal {
        Scal D = yc;
        Scal H = yh;
        H /= D;
        using std::sqrt;
        auto sech = [](Scal t) { return 1. / std::cosh(t); };

        Scal x = xx[0] - xc;
        Scal y = xx[1];
        x /= D;
        y /= D;
        Scal z = x * sqrt(3 * H / 4) * (1 - 5 * H / 8);
        Scal s = sech(z) * sech(z);
        Scal h = 1 + H * s - 4 * H * H * s * (1 - s) / 3;
        h *= D;

        return xx[1] - h;
      };
      for (auto c : m.Cells()) {
        auto x = m.GetCenter(c);
        Vect h = m.GetCellSize();
        fc[c] = GetLevelSetVolume<Scal, M::dim>(f, x, h);
      }
    };
  } else if (v == "solitoncos") {
    Scal xc(par.Double["soliton_xc"]);
    Scal yc(par.Double["soliton_yc"]);
    Scal yh(par.Double["soliton_yh"]);
    Scal xw(par.Double["soliton_xw"]);
    return [xc, yc, yh, xw](FieldCell<Scal>& fc, const M& m) {
      auto f = [xc, yh, yc, xw](const Vect& xx) -> Scal {
        Scal x = xx[0] - xc;
        Scal L = xw;
        Scal D = yc;
        Scal p2 = 2 * M_PI;
        Scal H = yh;
        Scal h = D + std::cos(p2 * x / L) * H;
        return xx[1] - h;
      };
      for (auto c : m.Cells()) {
        auto x = m.GetCenter(c);
        Vect h = m.GetCellSize();
        fc[c] = GetLevelSetVolume<Scal, M::dim>(f, x, h);
      }
    };
  } else if (v == "wavelamb") {
    Scal a0(par.Double["wavelamb_a0"]);
    Scal xc(par.Double["wavelamb_xc"]);
    Scal h(par.Double["wavelamb_h"]);
    Scal k(par.Double["wavelamb_k"]);

    return [a0, xc, h, k](FieldCell<Scal>& fc, const M& m) {
      auto f = [a0, xc, h, k](const Vect& xx) -> Scal {
        using std::tanh;
        using std::sin;
        using std::cos;
        using std::pow;
        Scal a = a0;
        Scal x = xx[0] - xc;
        Scal y = xx[1] - h;

        Scal eps = a * k;
        Scal chi = 1.0 / tanh(h * k);
        Scal eta = (1.0 / 4.0) * a * chi * eps * (3 * pow(chi, 2) - 1) *
                       cos(2 * k * x) +
                   a * pow(eps, 2) *
                       ((1.0 / 64.0) *
                            (24 * pow(chi, 6) + 3 * pow(pow(chi, 2) - 1, 2)) *
                            cos(3 * k * x) +
                        (1.0 / 8.0) * (-3 * pow(chi, 4) + 9 * pow(chi, 2) - 9) *
                            cos(k * x)) +
                   a * cos(k * x);

        return y - eta;
      };
      for (auto c : m.Cells()) {
        fc[c] =
            GetLevelSetVolume<Scal, M::dim>(f, m.GetCenter(c), m.GetCellSize());
      }
    };
  } else if (v == "solitonwang") {
    Scal xc(par.Double["soliton_xc"]);
    Scal yc(par.Double["soliton_yc"]);
    Scal yh(par.Double["soliton_yh"]);
    Scal xw(par.Double["soliton_xw"]);
    Scal e(par.Double["soliton_eps"]);
    return [xc, yc, yh, xw, e](FieldCell<Scal>& fc, const M& m) {
      auto f = [xc, yh, yc, xw, e](const Vect& xx) -> Scal {
        using std::cos;
        Scal a = yh;
        Scal la = xw;
        Scal k = 2. * M_PI / la;
        Scal x = xx[0] - xc;
        Scal kx = k * x;
        Scal h = yc + a / la *
                          (cos(kx) + 0.5 * e * cos(2 * kx) +
                           3. / 8 * sqr(e) * cos(3. * kx));
        return xx[1] - h;
      };
      for (auto c : m.Cells()) {
        auto x = m.GetCenter(c);
        Vect h = m.GetCellSize();
        fc[c] = GetLevelSetVolume<Scal, M::dim>(f, x, h);
      }
    };
  } else if (v == "box") {
    Vect xc(par.Vect["box_c"]);
    Scal s = par.Double["box_s"];
    return [xc, s](FieldCell<Scal>& fc, const M& m) {
      for (auto c : m.Cells()) {
        fc[c] = (xc - m.GetCenter(c)).norminf() < s * 0.5 ? 1. : 0.;
      }
    };
  } else if (v == "line") {
    Vect xc(par.Vect["line_c"]); // center
    Vect n(par.Vect["line_n"]); // normal
    Scal h(par.Double["line_h"]); // thickness
    return [xc, n, h](FieldCell<Scal>& fc, const M& m) {
      for (auto c : m.Cells()) {
        Scal d = (m.GetCenter(c) - xc).dot(n);
        fc[c] = 1. / (1. + std::exp(-d / h));
      }
    };
  } else if (v == "sin") {
    Vect k;
    if (auto p = par.Vect.Find("sin_k")) {
      k = Vect(*p);
    } else {
      k = Vect(2. * M_PI);
    }

    return [k](FieldCell<Scal>& fc, const M& m) {
      auto sin = [](Vect x) {
        Vect res;
        for (auto d : M::dirs) {
          res[d] = std::sin(x[d]);
        }
        return res;
      };
      for (auto c : m.Cells()) {
        fc[c] = sin(m.GetCenter(c) * k).prod();
      }
    };
  } else if (v == "sinc") {
    Vect k(par.Vect["sinc_k"]);
    return [k](FieldCell<Scal>& fc, const M& m) {
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
  } else if (v == "grid") { // see init_cl.h for grid of different colors
    return [](FieldCell<Scal>& fc, const M& m) {
      for (auto c : m.Cells()) {
        fc[c] = 1;
      }
    };
  } else if (v == "readplain") {
    return [](FieldCell<Scal>& fc, const M& m) { fc.Reinit(m, 1); };
  } else if (v == "zero") {
    return [](FieldCell<Scal>& fc, const M& m) {
      for (auto c : m.Cells()) {
        fc[c] = 0;
      }
    };
  } else {
    fassert(false, "Unknown init_vf=" + v);
  }
  return std::function<void(FieldCell<Scal>&, const M&)>();
}

template <class M>
void InitVf(
    FieldCell<typename M::Scal>& fcu, const Vars& var, M& m, bool verbose) {
  auto sem = m.GetSem("initvf");
  struct {
    std::vector<char> buf;
  } * ctx(sem);

  if (sem("zero")) {
    fcu.Reinit(m, 0);
  }
  std::string v = var.String["init_vf"];
  if (v == "list") {
    if (sem("list-bcast")) {
      if (m.IsRoot()) {
        std::stringstream path(var.String["list_path"]);
        std::string fname;
        path >> fname;
        if (fname == "inline") {
          if (verbose) {
            std::cerr
                << "InitVf: Reading inline list of primitives from list_path"
                << std::endl;
          }
          ctx->buf = std::vector<char>(
              std::istreambuf_iterator<char>(path),
              std::istreambuf_iterator<char>());
        } else {
          std::ifstream fin(fname);
          if (verbose) {
            std::cerr << "InitVf: Open list of primitives '" << fname << "'"
                      << std::endl;
          }
          fassert(fin.good(), "Can't open list of primitives '" + fname + "'");
          ctx->buf = std::vector<char>(
              std::istreambuf_iterator<char>(fin),
              std::istreambuf_iterator<char>());
        }
      }
      m.Bcast(&ctx->buf);
    }
    if (sem("list-local")) {
      std::stringstream list;
      std::copy(
          ctx->buf.begin(), ctx->buf.end(), std::ostream_iterator<char>(list));
      InitVfList(fcu, list, var.Int["list_ls"], var.Int["dim"], m, verbose);
    }
  } else if (v == "hdf") {
    if (sem.Nested()) {
      Hdf<M>::Read(fcu, var.String["init_vf_hdf_path"], m);
    }
  } else {
    if (sem("local")) {
      auto func = CreateInitU<M>(var, m.IsRoot());
      func(fcu, m);
    }
  }
  if (sem("comm")) {
    m.Comm(&fcu);
  }
}
