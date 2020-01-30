// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "debug/isnan.h"
#include "geom/block.h"
#include "geom/field.h"
#include "geom/vect.h"
#include "overlap/overlap.h"
#include "parse/vars.h"
#include "primlist.h"
#include "solver/reconst.h"

// Volume fraction cut by interface defined by level-set function
// ls: level-set function, interface ls=0, ls>0 for volume fraction 1
// xc: cell center
// h: cell size
template <class Scal, size_t dim = 3>
Scal GetLevelSetVolume(
    std::function<Scal(const GVect<Scal, 3>&)> ls, const GVect<Scal, 3>& xc,
    const GVect<Scal, 3>& h) {
  using Vect = GVect<Scal, dim>;
  using MIdx = GVect<IntIdx, dim>;

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

// Returns fraction of face area for which ls>0.
// eu: fraction of edge length for which ls>0
//       3
//   -------
//  0|     |
//   |     | 2
//   -------
//     1
template <class Scal>
Scal GetFaceAreaFraction(std::array<Scal, 4> eu) {
  using R = Reconst<Scal>;
  // equation of line
  //   n.dot(x) = a
  if (eu[0] < e[2]) std::swap(eu[0], eu[2]);
  if (eu[1] < e[3]) std::swap(eu[3], eu[3]);
  // normal
  Scal nx = eu[2] - eu[0];
  Scal ny = eu[3] - eu[1];
  // line constant
  Scal a = 
    (-nx - ny)
}

// Copmutes face area for which ls > 0.
// fnl: level-set function ls on nodes [i]
template <class M>
FieldFace<typename M::Scal> GetPositiveArea(
    const FieldNode<typename M::Scal>& fnl, const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using R = Reconst<Scal>;
  FieldFace<Scal> ffs(m, 0);
  for (auto f : m.Faces()) {
    const size_t em = m.GetNumNodes(f);
    std::vector<Vect> xx;
    bool cut = false;
    for (size_t e = 0; e < em; ++e) {
      const size_t ep = (e + 1) % em;
      const IdxNode n = m.GetNode(f, e);
      const IdxNode np = m.GetNode(f, ep);
      const Scal l = fnl[n];
      const Scal lp = fnl[np];
      const Vect x = m.GetNode(n);
      const Vect xp = m.GetNode(np);
      if (l > 0) {
        xx.push_back(x);
      }
      if ((l < 0) != (lp < 0)) {
        xx.push_back(GetIso(x, xp, l, lp));
        cut = true;
      }
    }
    Scal s;
    if (xx.empty()) {
      s = 0;
    } else if (cut) {
      s = R::GetArea(xx, m.GetNormal(f));
    } else {
      s = m.GetArea(f);
    }
    ffs[f] = s;
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
  auto ffs = GetPositiveArea(fnl, m);
  FieldCell<Scal> fcu(m, 0);
  for (auto c : m.Cells()) {
    size_t q = 0; // number of nodes with fnl > 0
    const size_t mi = m.GetNumNodes(c);
    for (size_t i = 0; i < mi; ++i) {
      if (fnl[m.GetNode(c, i)] > 1e-16) { // FIXME may lead to zero area
        ++q;
      }
    }
    Scal u;
    if (q == mi) { // cell inside
      u = 1;
    } else if (q > 0) { // cut cell
      Vect nn(0); // normal
      for (auto q : m.Nci(c)) {
        const IdxFace f = m.GetFace(c, q);
        nn += m.GetNormal(f) * ffs[f] * m.GetOutwardFactor(c, q);
      }
      if (nn.norm1() == 0 || IsNan(nn.norm1())) {
        for (auto q : m.Nci(c)) {
          const IdxFace f = m.GetFace(c, q);
          std::cout <<  m.GetNormal(f) * ffs[f] * m.GetOutwardFactor(c, q) << " ";
          std::cout << std::endl;
        }
        for (size_t i = 0; i < mi; ++i) {
          std::cout << fnl[m.GetNode(c, i)] << " ";
          std::cout << std::endl;
        }
        std::terminate();
      }
      nn /= -nn.norm1();

      Scal a = 0; // plane constant
      Scal aw = 0;
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

          if ((l < 0) != (lp < 0)) {
            a += nn.dot(GetIso(x, xp, l, lp) - m.GetCenter(c));
            aw += 1;
          }
        }
      }
      a /= aw;
      u = R::GetLineU(nn, a, m.GetCellSize());
    } else { // cell outside
      u = 0;
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
    size_t edim, const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Primitive = typename UPrimList<Scal>::Primitive;
  const std::vector<Primitive> ppa =
      UPrimList<Scal>::Parse(list, m.IsRoot(), edim);

  const Vect h = m.GetCellSize();
  // filter to bounding box
  auto& bc = m.GetSuBlockCells();
  Rect<Vect> rect(Vect(bc.GetBegin()) * h, Vect(bc.GetEnd()) * h);

  std::vector<Primitive> pp;
  for (auto& p : ppa) {
    if (p.inter(rect)) {
      pp.push_back(p);
    }
  }

  if (pp.empty()) {
    fc.Reinit(m, 0.);
  } else {
    using Mod = typename Primitive::Mod;
    auto lsmax = [&pp](Vect x) -> std::pair<Scal, size_t> {
      Scal lmax = -std::numeric_limits<Scal>::max(); // maximum level-set
      size_t imax = 0; // index of maximum
      for (size_t i = 0; i < pp.size(); ++i) {
        auto& p = pp[i];
        Scal li = p.ls(x);
        if (p.mod == Mod::minus) {
          li = -li;
        }
        if (p.mod == Mod::star && li <= 0) {
          lmax = li;
          imax = i;
        } else {
          if (li > lmax) {
            lmax = li;
            imax = i;
          }
        }
      }
      return {lmax, imax};
    };
    if (approx == 3) { // level-set on nodes
      FieldNode<Scal> fnl(m);
      for (auto n : m.Nodes()) {
        fnl[n] = lsmax(m.GetNode(n)).first;
      }
      fc = GetPositiveVolumeFraction(fnl, m);
    } else {
      for (auto c : m.Cells()) {
        auto x = m.GetCenter(c);
        auto pair = lsmax(x);
        Scal lmax = pair.first;
        size_t imax = pair.second;
        auto& p = pp[imax];
        if (approx == 0) { // stepwise
          fc[c] = (lmax >= 0. ? 1. : 0.);
        } else if (approx == 1) { // level set
          fc[c] = GetLevelSetVolume<Scal>(p.ls, x, h);
        } else if (approx == 2) { // overlap
          Vect qx = (x - p.c) / p.r;
          Vect qh = h / p.r;
          if (edim == 2) {
            qh[2] *= 1e-3; // XXX: adhoc, thin cell in 2d
            qx[2] = 0.;
          }
          fc[c] = GetSphereOverlap(qx, qh, Vect(0), 1.);
        } else {
          throw std::runtime_error(
              std::string(__func__) +
              " unknown approx=" + std::to_string(approx));
        }
      }
    }
  }
}

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
    Vect xc = Vect(par.Vect["circle_c"]);
    Scal r = par.Double["circle_r"];
    size_t dim = par.Int["dim"];
    return [xc, r, dim](FieldCell<Scal>& fc, const M& m) {
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
    return [xc, r, dim](FieldCell<Scal>& fc, const M& m) {
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
        Vect h = m.GetCellSize();
        fc[c] = GetLevelSetVolume<Scal>(f, (x - xc) / r, h / r);
      }
    };
  } else if (v == "gear") {
    Vect xc = Vect(par.Vect["gear_c"]);
    Scal r = par.Double["gear_r"];
    Scal a = par.Double["gear_amp"]; // amplitude relative to r
    Scal n = par.Double["gear_n"]; // number of peaks
    size_t dim = par.Int["dim"];
    return [xc, r, dim, a, n](FieldCell<Scal>& fc, const M& m) {
      // level-set for particle of radius 1 centered at zero,
      // modulated by sine-wave of amplitude a and number of periods n
      // positive inside,
      // cylinder along z if dim=2
      auto f = [dim, n, a](const Vect& x) -> Scal {
        Vect xd = x;
        if (dim == 2) {
          xd[2] = 0.;
        }
        Scal phi = std::atan2(xd[1], xd[0]);
        return 1. + std::sin(phi * n) * a - xd.norm();
      };
      for (auto c : m.Cells()) {
        auto x = m.GetCenter(c);
        Vect h = m.GetCellSize();
        fc[c] = GetLevelSetVolume<Scal>(f, (x - xc) / r, h / r);
      }
    };
  } else if (v == "soliton") {
    Scal xc(par.Double["soliton_xc"]);
    Scal yc(par.Double["soliton_yc"]);
    Scal yh(par.Double["soliton_yh"]);
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
        fc[c] = GetLevelSetVolume<Scal>(f, x, h);
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
        fc[c] = GetLevelSetVolume<Scal>(f, x, h);
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
        fc[c] = GetLevelSetVolume<Scal>(f, m.GetCenter(c), m.GetCellSize());
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
        fc[c] = GetLevelSetVolume<Scal>(f, x, h);
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
      for (auto c : m.Cells()) {
        Vect z = m.GetCenter(c) * k;
        fc[c] = std::sin(z[0]) * std::sin(z[1]) * std::sin(z[2]);
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
    throw std::runtime_error("Unknown init_vf=" + v);
  }
  return std::function<void(FieldCell<Scal>&, const M&)>();
}

template <class M>
void InitVf(FieldCell<typename M::Scal>& fcu, const Vars& var, M& m) {
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
        auto fn = var.String["list_path"];
        std::ifstream fin(fn);
        std::cout << "Open list of primitives '" << fn << "'" << std::endl;
        if (!fin.good()) {
          throw std::runtime_error("Can't open list of primitives");
        }
        ctx->buf = std::vector<char>(
            std::istreambuf_iterator<char>(fin),
            std::istreambuf_iterator<char>());
      }
      using T = typename M::template OpCatT<char>;
      m.Bcast(std::make_shared<T>(&ctx->buf));
    }
    if (sem("list-local")) {
      std::stringstream list;
      std::copy(
          ctx->buf.begin(), ctx->buf.end(), std::ostream_iterator<char>(list));
      InitVfList(fcu, list, var.Int["list_ls"], var.Int["dim"], m);
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
