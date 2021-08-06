// Created by Petr Karnakov on 28.12.2019
// Copyright 2019 ETH Zurich

#include "init.h"
#include "geom/mesh.h"
#include "init_cl.h"
#include "init_sig.h"
#include "init_u.h"
#include "overlap/overlap.h"

template <class M>
void InitVfList(
    FieldCell<typename M::Scal>& fc, std::istream& primlist, int approx,
    size_t edim, const M& m, bool verbose) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Primitive = typename UPrimList<Vect>::Primitive;
  const std::vector<Primitive> ppa =
      UPrimList<Vect>::GetPrimitives(primlist, edim);
  if (verbose && m.IsRoot()) {
    std::cerr << "Read " << ppa.size() << " primitives" << std::endl;
  }

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
    auto lsmax0 = [&pp](Vect x) -> std::pair<Scal, size_t> {
      Scal lmax = -std::numeric_limits<Scal>::max(); // maximum level-set
      size_t imax0 = 0; // index of maximum
      for (size_t i = 0; i < pp.size(); ++i) {
        auto& p = pp[i];
        Scal li = p.ls(x);
        if (p.mod_minus) {
          li = -li;
        }
        if (p.mod_and) {
          lmax = std::min(lmax, li);
        } else {
          if (li > lmax) {
            lmax = li;
            imax0 = i;
          }
        }
      }
      return {lmax, imax0};
    };
    if (approx == 0) { // stepwise
      for (auto c : m.Cells()) {
        fc[c] = (lsmax0(m.GetCenter(c)).first >= 0 ? 1 : 0);
      }
    } else if (approx == 1) { // level set
      for (auto c : m.Cells()) {
        const auto x = m.GetCenter(c);
        const auto& p = pp[lsmax0(x).second];
        fc[c] = GetLevelSetVolume<Scal>(p.ls, x, h);
      }
    } else if (approx == 2) { // overlap
#if USEFLAG(OVERLAP)
      for (auto c : m.Cells()) {
        const auto x = m.GetCenter(c);
        const auto& p = pp[lsmax0(m.GetCenter(c)).second];
        Vect qx = (x - p.c) / p.r;
        Vect qh = h / p.r;
        if (edim == 2 && M::dim == 3) {
          qx[2] = 0.;
          qh[2] *= 1e-3; // XXX: adhoc, thin cell for 2d
        }
        using Vect3 = generic::Vect<Scal, 3>;
        fc[c] = GetSphereOverlap(Vect3(qx), Vect3(qh), Vect3(0), 1);
      }
#else
      fassert(false, "overlap is disabled");
#endif
    } else if (approx == 3) { // level-set on nodes
      FieldNode<Scal> fnl(m);
      for (auto n : m.Nodes()) {
        fnl[n] = lsmax0(m.GetNode(n)).first;
      }
      fc = GetPositiveVolumeFraction(fnl, m);
    } else {
      fassert(false, "unknown approx=" + std::to_string(approx));
    }
  }
}

template <class M>
void InitOverlappingComponents(
    std::istream& primlist, const Multi<FieldCell<typename M::Scal>*>& fcu,
    const Multi<FieldCell<typename M::Scal>*>& fccl,
    const GRange<size_t>& layers, const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Primitive = typename UPrimList<Vect>::Primitive;
  constexpr Scal kClNone = -1; // no color
  const std::vector<Primitive> primitives =
      UPrimList<Vect>::GetPrimitives(primlist, m.flags.edim);
  fassert_equal(fcu.size(), layers.size());
  fassert_equal(fccl.size(), layers.size());
  for (auto l : layers) {
    fcu[l]->Reinit(m, 0);
    fccl[l]->Reinit(m, kClNone);
  }
  const auto h = m.GetCellSize();
  auto levelset = [&](size_t i, Vect x) {
    auto& p = primitives[i];
    return p.ls(x) * (p.mod_minus ? -1 : 1);
  };
  for (auto c : m.AllCellsM()) {
    const size_t none = -1;
    size_t imax0 = none; // primitive with largest levelset in c
    size_t imax1 = none; // primitive with second largest levelset in c
    Scal lsmax0 = -std::numeric_limits<Scal>::max();
    Scal lsmax1 = -std::numeric_limits<Scal>::max();
    Vect image0(0); // image vector for periodic conditions
    Vect image1(0); // image vector for periodic conditions

    const auto wimages = [&]() { //
      using MIdx = generic::MIdx<M::dim>;
      MIdx begin(0);
      MIdx size(1);
      for (auto d : m.dirs) {
        if (m.flags.is_periodic[d]) {
          begin[d] = -1;
          size[d] = 3;
        }
      }
      return GBlock<int, M::dim>{begin, size};
    }();

    for (size_t i = 0; i < primitives.size(); ++i) {
      for (auto wimage : wimages) {
        const Vect image = Vect(wimage) * m.GetGlobalLength();
        const Scal ls = levelset(i, c.center() + image);
        if (ls > lsmax0) {
          imax0 = i;
          lsmax0 = ls;
          image0 = image;
        }
      }
    }
    for (size_t i = 0; i < primitives.size(); ++i) {
      for (auto wimage : wimages) {
        const Vect image = Vect(wimage) * m.GetGlobalLength();
        const Scal ls = levelset(i, c.center() + image);
        if (ls > lsmax1 && i != imax0) {
          imax1 = i;
          lsmax1 = ls;
          image1 = image;
        }
      }
    }
    fassert(imax0 != none);
    const Scal umax0 =
        (imax0 == none ? 0
                       : GetLevelSetVolume<Scal, M::dim>(
                             primitives[imax0].ls, c.center() + image0, h));
    const Scal umax1 =
        (imax1 == none ? 0
                       : GetLevelSetVolume<Scal, M::dim>(
                             primitives[imax1].ls, c.center() + image1, h));
    if (imax1 == none) {
      (*fccl[0])[c] = imax0;
      (*fcu[0])[c] = umax0;
    } else if (umax0 + umax1 >= 1) { // two or more positive level-sets
      auto f0 = [&](const Vect& x) -> Scal {
        return levelset(imax0, x + image0) - levelset(imax1, x + image1);
      };
      auto f1 = [&](const Vect& x) -> Scal {
        return levelset(imax1, x + image1) - levelset(imax0, x + image0);
      };
      const Scal u0 = GetLevelSetVolume<Scal, M::dim>(f0, c.center, h);
      const Scal u1 = GetLevelSetVolume<Scal, M::dim>(f1, c.center, h);
      if (u0 > 0) {
        for (auto l : layers) {
          if ((*fccl[l])[c] == kClNone) {
            (*fccl[l])[c] = imax0;
            (*fcu[l])[c] = u0;
            break;
          }
        }
      }
      if (u1 > 0) {
        for (auto l : layers) {
          if ((*fccl[l])[c] == kClNone) {
            (*fccl[l])[c] = imax1;
            (*fcu[l])[c] = u1;
            break;
          }
        }
      }
    } else if (umax0 > 0) {
      (*fccl[0])[c] = imax0;
      (*fcu[0])[c] = umax0;
    } else if (umax1 > 0) {
      (*fccl[0])[c] = imax1;
      (*fcu[0])[c] = umax1;
    }
    // Override with primitives having `mod_and == true`
    for (size_t i = 0; i < primitives.size(); ++i) {
      if (primitives[i].mod_and) {
        for (auto wimage : wimages) {
          const Vect image = Vect(wimage) * m.GetGlobalLength();
          const Scal u = GetLevelSetVolume<Scal, M::dim>(
              primitives[i].ls, c.center() + image, h);
          if (u > 0) {
            if (u == 1) {
              for (auto l : layers) {
                (*fccl[l])[c] = kClNone;
                (*fcu[l])[c] = 0;
              }
            }
            Scal sum = 0; // current sum of existing volume fractions
            for (auto l : layers) {
              if ((*fccl[l])[c] != kClNone) {
                sum += (*fcu[l])[c];
              }
            }
            // Reduce existing volume fractions to make space for the new one
            if (sum > 0) {
              for (auto l : layers) {
                if ((*fccl[l])[c] != kClNone) {
                  (*fcu[l])[c] *= (1 - u) / sum;
                }
              }
            }
            // Add new component
            for (auto l : layers) {
              if ((*fccl[l])[c] == kClNone) {
                (*fccl[l])[c] = i;
                (*fcu[l])[c] = u;
                break;
              }
            }
          }
        }
      }
      imax1 = none;
    }
  }
}

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

std::stringstream ReadPrimList(std::string list_path, bool verbose) {
  std::stringstream path(list_path);
  std::stringstream res;
  std::string fname;
  path >> fname;
  if (fname == "inline") {
    if (verbose) {
      std::cerr << "Reading inline list of primitives from list_path"
                << std::endl;
    }
    res << path.rdbuf();
  } else {
    std::ifstream fin(fname);
    if (verbose) {
      std::cerr << "Open list of primitives '" << fname << "'" << std::endl;
    }
    fassert(fin.good(), "Can't open list of primitives '" + fname + "'");
    res << fin.rdbuf();
  }
  return res;
}

#define XX(M)                                                                 \
  template void InitVf(                                                       \
      FieldCell<typename M::Scal>& fcu, const Vars& var, M& m, bool verbose); \
  template std::function<void(                                                \
      FieldCell<typename M::Scal>&, const FieldCell<typename M::Scal>&,       \
      const M&)>                                                              \
  CreateInitCl(const Vars& par, bool verb);                                   \
  template std::function<void(FieldCell<typename M::Scal>&, const M&)>        \
  CreateInitU(const Vars& par, bool verb);                                    \
  template std::function<void(FieldCell<typename M::Scal>&, const M&)>        \
  CreateInitSig(const Vars& var);                                             \
  template void InitOverlappingComponents(                                    \
      std::istream& primlist, const Multi<FieldCell<typename M::Scal>*>& fcu, \
      const Multi<FieldCell<typename M::Scal>*>& fccl,                        \
      const GRange<size_t>& layers, const M& m);
#define COMMA ,
#define X(dim) XX(MeshCartesian<double COMMA dim>)
MULTIDIMX
