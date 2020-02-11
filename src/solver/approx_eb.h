// Created by Petr Karnakov on 11.02.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <array>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "cond.h"
#include "embed.h"
#include "geom/mesh.h"
#include "solver.h"

template <class Scal>
struct ULinear {
  static constexpr size_t dim = 3;
  using Vect = generic::Vect<Scal, dim>;

  template <size_t N>
  static std::array<Scal, N> Mul(
      std::array<Scal, N * N> a, std::array<Scal, N> x) {
    using Int = size_t;
    std::array<Scal, N> r;
    for (Int i = 0; i < N; ++i) {
      r[i] = 0;
      for (Int j = 0; j < N; ++j) {
        r[i] += a[i * N + j] * x[j];
      }
    }
    return r;
  }

  // Solves linear system a*x=b.
  template <size_t N>
  static std::array<Scal, N> SolveLinear(
      std::array<Scal, N * N> a, std::array<Scal, N> b) {
    using Int = size_t;
    auto aa = [&a](Int i, Int j) -> Scal& { return a[i * N + j]; };
    auto swaprows = [&aa, &b](Int i, Int ip) {
      if (i == ip) {
        return;
      }
      for (Int j = 0; j < N; ++j) {
        std::swap(aa(i, j), aa(ip, j));
      }
      std::swap(b[i], b[ip]);
    };
    auto ipivot = [&aa](const Int j) {
      Int imax = j;
      for (Int i = j + 1; i < N; ++i) {
        if (std::abs(aa(i, j)) > std::abs(aa(imax, j))) {
          imax = i;
        }
      }
      return imax;
    };
    auto addrow = [&aa, &b](Int i, Int ip, Scal ap) {
      for (Int j = 0; j < N; ++j) {
        aa(i, j) += aa(ip, j) * ap;
      }
      b[i] += b[ip] * ap;
    };
    for (Int j = 0; j < N; ++j) {
      const Int ip = ipivot(j);
      swaprows(ip, j);
      for (Int i = j + 1; i < N; ++i) {
        addrow(i, j, -aa(i, j) / aa(j, j));
      }
    }
    std::array<Scal, N> x;
    for (Int i = N; i > 0;) {
      --i;
      Scal t = b[i];
      for (Int j = i + 1; j < N; ++j) {
        t -= aa(i, j) * x[j];
      }
      x[i] = t / aa(i, i);
    }
    return x;
  }

  // Fits linear function to set of points and values
  //   u = g.dot(x) + u0
  // Returns {g, u0}.
  static std::pair<Vect, Scal> FitLinear(
      const std::vector<Vect>& xx, const std::vector<Scal>& uu) {
    assert(xx.size() == uu.size());
    // sum 0.5 * [ (g.dot(x[k]) + u0 - u[k]) ** 2 ] -> min
    using Int = size_t;
    constexpr Int N = dim + 1;
    std::array<Scal, N * N> a;
    std::array<Scal, N> b;
    auto aa = [&a](Int i, Int j) -> Scal& { return a[i * N + j]; };
    std::fill(a.begin(), a.end(), 0);
    std::fill(b.begin(), b.end(), 0);
    for (size_t k = 0; k < xx.size(); ++k) {
      for (Int i = 0; i < dim; ++i) {
        for (Int j = 0; j < dim; ++j) {
          aa(i, j) += xx[k][j] * xx[k][i];
        }
        aa(i, dim) += xx[k][i];
        b[i] += uu[k] * xx[k][i];
      }
      for (Int j = 0; j < dim; ++j) {
        aa(dim, j) += xx[k][j];
      }
      aa(dim, dim) += 1;
      b[dim] += uu[k];
    }

    auto v = SolveLinear(a, b);
    return {Vect(v[0], v[1], v[2]), v[3]};
  }

  // Fits linear function to set of points and values
  //   u = g.dot(x) + u0
  // Returns {g, u0}.
  static std::pair<generic::Vect<Vect, dim>, Vect> FitLinear(
      const std::vector<Vect>& xx, const std::vector<Vect>& uu) {
    std::pair<generic::Vect<Vect, dim>, Vect> p;
    for (size_t d = 0; d < dim; ++d) {
      std::vector<Scal> uud;
      for (auto u : uu) {
        uud.push_back(u[d]);
      }
      auto pd = FitLinear(xx, uud);
      p.second[d] = pd.second;
      p.first[0][d] = pd.first[0];
      p.first[1][d] = pd.first[1];
      p.first[2][d] = pd.first[2];
    }
    return p;
  }

  template <class T>
  static T EvalLinear(
      const std::pair<generic::Vect<T, dim>, T>& p,
      const generic::Vect<Scal, dim>& x) {
    auto& g = p.first;
    auto& u0 = p.second;
    return g[0] * x[0] + g[1] * x[1] + g[2] * x[2] + u0;
  }
};

template <class M_>
struct UEmbed {
  using M = M_;
  using EB = Embed<M>;
  static constexpr size_t dim = M::dim;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Type = typename EB::Type;

  struct CondEmbed {
    enum class Type { value, gradient };
    Scal u; // value or normal gradient
  };

  static FieldNode<Scal> InitEmbed(const M& m, const Vars& var, bool verb) {
    FieldNode<Scal> fnl(m); // level-set
    const auto name = var.String["eb_init"];
    if (name == "none") {
      fnl.Reinit(m, 1);
    } else if (name == "box") {
      const Vect xc(var.Vect["eb_box_c"]);
      const Vect r(var.Vect["eb_box_r"]);
      const Scal angle = M_PI * var.Double["eb_box_angle"];
      for (auto n : m.AllNodes()) {
        const Vect x = m.GetNode(n);
        auto rot = [angle](Vect xx) {
          const Scal sin = std::sin(angle);
          const Scal cos = std::cos(angle);
          const Scal x = xx[0];
          const Scal y = xx[1];
          const Scal z = xx[2];
          return Vect(x * cos - y * sin, x * sin + y * cos, z);
        };
        fnl[n] =
            (1 - (rot(x - xc) / r).norminf()) * (r / m.GetCellSize()).min();
      }
    } else if (name == "sphere") {
      const Vect xc(var.Vect["eb_sphere_c"]);
      const Vect r(var.Vect["eb_sphere_r"]);
      const Scal angle = M_PI * var.Double["eb_sphere_angle"];
      for (auto n : m.AllNodes()) {
        const Vect x = m.GetNode(n);
        auto rot = [angle](Vect xx) {
          const Scal sin = std::sin(angle);
          const Scal cos = std::cos(angle);
          const Scal x = xx[0];
          const Scal y = xx[1];
          const Scal z = xx[2];
          return Vect(x * cos - y * sin, x * sin + y * cos, z);
        };
        fnl[n] = (rot(x - xc) / r).norm() - 1;
      }
    } else if (name == "list") {
      // TODO revise with bcast
      const std::string fn = var.String["eb_list_path"];
      const size_t edim = var.Int["dim"];
      std::ifstream fin(fn);
      if (verb) {
        std::cout << "Open list of primitives '" << fn << "' for embed"
                  << std::endl;
      }
      if (!fin.good()) {
        throw std::runtime_error("Can't open list of primitives");
      }
      auto pp = UPrimList<Scal>::Parse(fin, verb, edim);

      for (auto n : m.AllNodes()) {
        Scal lmax = -std::numeric_limits<Scal>::max(); // maximum level-set
        for (auto& p : pp) {
          lmax = std::max(lmax, p.ls(m.GetNode(n)));
        }
        fnl[n] = lmax;
      }
    } else {
      throw std::runtime_error("Unknown eb_init=" + name);
    }

    if (var.Int["eb_init_inverse"]) {
      for (auto n : m.AllNodes()) {
        fnl[n] = -fnl[n];
      }
    }
    return fnl;
  }

  // feu: field on embedded boundaries [a]
  // Returns:
  // field on cells [a]
  template <class T>
  static FieldCell<T> Interpolate(const FieldEmbed<T>& feu, const EB& eb) {
    auto& m = eb.GetMesh();
    FieldCell<T> fcu(eb, T(0));
    for (auto c : eb.AllCells()) {
      if (eb.GetType(c) == Type::regular) {
        T sum(0);
        Scal sumw = 0;
        for (auto q : m.Nci(c)) {
          IdxFace f = eb.GetFace(c, q);
          sum += feu[f];
          sumw += 1.;
        }
        fcu[c] = sum / sumw;
      } else { // Type::cut
        const Scal w = 1 / std::abs(GetFaceOffset(c));
        T sum = feu[c] * w;
        Scal sumw = w;
        for (auto q : eb.Nci(c)) {
          IdxFace f = eb.GetFace(c, q);
          const Scal w = 1 / std::abs(GetFaceOffset(c, q));
          sum += feu[f] * w;
          sumw += w;
        }
        fcu[c] = sum / sumw;
      }
    }
    return fcu;
  }
  // Mid-point interpolation.
  // fcu: field [a]
  // bc: boundary conditions type, 0: value, 1: grad
  // bcv: value or normal gradient (grad dot GetNormal)
  // Returns:
  // field on embedded boundaries [s]
  template <class T>
  static FieldEmbed<T> Interpolate(
      const FieldCell<T>& fcu, const MapCondFace& mfc, size_t bc, T bcv,
      const EB& eb) {
    FieldEmbed<T> feu(eb, T(0));
    for (auto f : eb.SuFaces()) {
      const IdxCell cm = eb.GetCell(f, 0);
      const IdxCell cp = eb.GetCell(f, 1);
      const Scal a = 0.5;
      feu[f] = fcu[cm] * (1 - a) + fcu[cp] * a;
    }
    for (auto c : eb.SuCFaces()) {
      if (bc == 0) {
        feu[c] = bcv;
      } else if (bc == 1) {
        feu[c] = fcu[c] + bcv * eb.GetFaceOffset(c);
      } else {
        throw std::runtime_error(
            "Interpolate: unknown bc=" + std::to_string(bc));
      }
    }
    InterpolateB(fcu, mfc, feu.GetFieldFace(), eb.GetMesh());
    return feu;
  }
  // Interpolation using gradient with iterative correction.
  // High order gradient-based interpolation with extrapolation
  // boundary conditions.
  // fcu: field [a]
  // bc: boundary conditions type, 0: value, 1: grad
  // bcv: value or normal gradient (grad dot GetNormal)
  // Returns:
  // field on embedded boundaries [s]
  template <class T>
  static FieldEmbed<T> InterpolateIterative(
      const FieldCell<T>& fcu, const EB& eb) {
    // interpolate from cells to faces with zero-derivative conditions
    auto feu = Interpolate(fcu, MapCondFace(), 1, T(0), eb);

    FieldCell<Vect> fcg(eb, Vect(0));

    for (size_t i = 0; i < 20; ++i) {
      // compute gradient from current interpolant
      fcg = GradientLinearFit(feu, eb);

      // goal: compute first-order accurate gradient
      // and using the gradient, do second-order accurate interpolation

      // update feu from fcu and fcg
      for (auto f : eb.SuFaces()) {
        const IdxCell cm = eb.GetCell(f, 0);
        const IdxCell cp = eb.GetCell(f, 1);
        if (eb.GetType(cm) == Type::regular &&
            eb.GetType(cp) == Type::regular) {
          feu[f] = (fcu[cm] + fcu[cp]) * 0.5;
        } else {
          const Vect xm = eb.GetCellCenter(cm);
          const Vect xp = eb.GetCellCenter(cp);
          const Vect xf = eb.GetFaceCenter(f);
          feu[f] = (fcu[cm] + fcg[cm].dot(xf - xm) + fcu[cp] +
                    fcg[cp].dot(xf - xp)) *
                   0.5;
        }
      }
      for (auto c : eb.SuCFaces()) {
        const Vect xc = eb.GetCellCenter(c);
        const Vect xf = eb.GetFaceCenter(c);
        feu[c] = fcu[c] + fcg[c].dot(xf - xc);
      }
    }
    return feu;
  }
  // Upwind interpolation.
  // fcu: field [a]
  // ffv: volume flux
  // mfc: conditions on faces
  // bc: boundary conditions type, 0: value, 1: grad
  // bcv: value or normal gradient (grad dot GetNormal)
  // Returns:
  // field on embedded boundaries [s]
  template <class T>
  static FieldEmbed<T> InterpolateUpwind(
      const FieldCell<T>& fcu, const FieldEmbed<Scal>& fev,
      const MapCondFace& mfc, size_t bc, T bcv, const EB& eb) {
    FieldEmbed<T> feu(eb, T(0));
    for (auto f : eb.SuFaces()) {
      const IdxCell c = (fev[f] > 0 ? eb.GetCell(f, 0) : eb.GetCell(f, 1));
      feu[f] = fcu[c];
    }
    for (auto c : eb.SuCFaces()) {
      if (fev[c] > 0) {
        feu[c] = fcu[c];
      } else {
        if (bc == 0) {
          feu[c] = bcv;
        } else if (bc == 1) {
          feu[c] = fcu[c] + bcv * eb.GetFaceOffset(c);
        } else {
          throw std::runtime_error(
              "Interpolate: unknown bc=" + std::to_string(bc));
        }
      }
    }
    InterpolateB(fcu, mfc, feu.GetFieldFace(), eb.GetMesh());
    return feu;
  }
  // Interpolates field from cells to faces with an upwind scheme,
  // bilinear in cut faces, and first order upwind in embed faces.
  // fcu: field [s]
  // fcg: gradient [s]
  // ffv: volume flux [i]
  // sc: interpolation scheme
  // Returns:
  // field on faces boundaries [s]
  static FieldFace<Scal> InterpolateUpwindBilinear(
      const FieldCell<Scal>& fcu, const FieldCell<typename M::Vect>& fcg,
      const FieldFace<Scal>& ffv, ConvSc sc, const EB& eb) {
    auto& m = eb.GetMesh();
    FieldFace<Scal> ffu(eb, 0.);
    // f = fmm*a[0] + fm*a[1] + fp*a[2]
    const std::array<Scal, 3> a = GetCoeff<Scal>(sc);
    for (auto f : eb.SuFaces()) {
      const IdxCell cm = eb.GetCell(f, 0);
      const IdxCell cp = eb.GetCell(f, 1);
      if (ffv[f] > 0) {
        ffu[f] = 4. * a[0] * fcg[cm].dot(m.GetVectToCell(f, 0)) +
                 a[1] * fcu[cm] + (a[2] + a[0]) * fcu[cp];
      } else if (ffv[f] < 0) {
        ffu[f] = 4. * a[0] * fcg[cp].dot(m.GetVectToCell(f, 1)) +
                 a[1] * fcu[cp] + (a[2] + a[0]) * fcu[cm];
      } else {
        ffu[f] = (fcu[cm] + fcu[cp]) * 0.5;
      }
    }
    return InterpolateBilinearFaces(ffu, eb);
  }
  static FieldEmbed<Scal> InterpolateUpwindBilinear(
      const FieldCell<Scal>& fcu, const FieldCell<typename M::Vect>& fcg,
      const MapCondFace& mfc, size_t bc, const MapCell<Scal>& mcu,
      const FieldEmbed<Scal>& fev, ConvSc sc, const EB& eb) {
    FieldEmbed<Scal> feu(eb, 0);
    feu.GetFieldFace() =
        InterpolateUpwindBilinear(fcu, fcg, fev.GetFieldFace(), sc, eb);
    InterpolateUpwindEmbedFaces(fcu, bc, mcu, fev, feu, eb);
    InterpolateB(fcu, mfc, feu.GetFieldFace(), eb.GetMesh());
    return feu;
  }
  static FieldEmbed<Scal> InterpolateUpwindBilinear(
      const FieldCell<Scal>& fcu, const FieldCell<typename M::Vect>& fcg,
      const MapCondFace& mfc, size_t bc, Scal bcv, const FieldEmbed<Scal>& fev,
      ConvSc sc, const EB& eb) {
    MapCell<Scal> mcu;
    for (auto c : eb.SuCFaces()) {
      mcu[c] = bcv;
    }
    return InterpolateUpwindBilinear(fcu, fcg, mfc, bc, mcu, fev, sc, eb);
  }
  // Upwind interpolation.
  // fcu: field [s]
  // fcg: gradient [s]
  // mfc: conditions on faces
  // bc: boundary conditions type, 0: value, 1: grad
  // bcv: value or normal gradient (grad dot GetNormal)
  // ffv: volume flux [i]
  // sc: interpolation scheme
  // Returns:
  // field on embedded boundaries [s]
  // TODO: interpolation from upwind cell on embedded boundaries
  //       (see InterpolateUpwind() above)
  static FieldEmbed<Scal> InterpolateUpwind(
      const FieldCell<Scal>& fcu, const FieldCell<typename M::Vect>& fcg,
      const MapCondFace& mfc, size_t bc, const MapCell<Scal>& mcu,
      const FieldEmbed<Scal>& fev, ConvSc sc, const EB& eb) {
    auto& m = eb.GetMesh();
    FieldEmbed<Scal> feu(eb, 0.);
    // f = fmm*a[0] + fm*a[1] + fp*a[2]
    const std::array<Scal, 3> a = GetCoeff<Scal>(sc);
    for (auto f : eb.SuFaces()) {
      const IdxCell cm = eb.GetCell(f, 0);
      const IdxCell cp = eb.GetCell(f, 1);
      if (eb.GetType(cm) == Type::regular && eb.GetType(cp) == Type::regular) {
        if (fev[f] > 0) {
          feu[f] = 4. * a[0] * fcg[cm].dot(m.GetVectToCell(f, 0)) +
                   a[1] * fcu[cm] + (a[2] + a[0]) * fcu[cp];
        } else if (fev[f] < 0) {
          feu[f] = 4. * a[0] * fcg[cp].dot(m.GetVectToCell(f, 1)) +
                   a[1] * fcu[cp] + (a[2] + a[0]) * fcu[cm];
        } else {
          feu[f] = (fcu[cm] + fcu[cp]) * 0.5;
        }
      } else {
        const IdxCell c = (fev[f] > 0 ? cm : cp);
        feu[f] = fcu[c];
      }
    }
    InterpolateUpwindEmbedFaces(fcu, bc, mcu, fev, feu, eb);
    InterpolateB(fcu, mfc, feu.GetFieldFace(), eb.GetMesh());
    return feu;
  }
  static FieldEmbed<Scal> InterpolateUpwind(
      const FieldCell<Scal>& fcu, const FieldCell<typename M::Vect>& fcg,
      const MapCondFace& mfc, size_t bc, Scal bcv, const FieldEmbed<Scal>& fev,
      ConvSc sc, const EB& eb) {
    MapCell<Scal> mcu;
    for (auto c : eb.SuCFaces()) {
      mcu[c] = bcv;
    }
    return InterpolateUpwind(fcu, fcg, mfc, bc, mcu, fev, sc, eb);
  }
  // feu: field on embedded boundaries [a]
  // Returns:
  // gradient on cells [a]
  static FieldCell<Vect> Gradient(const FieldEmbed<Scal>& feu, const EB& eb) {
    auto& m = eb.GetMesh();
    FieldCell<Vect> fcg(eb, Vect(0));
    for (auto c : eb.AllCells()) {
      if (eb.GetType(c) == Type::regular) {
        Vect sum(0);
        for (auto q : m.Nci(c)) {
          const IdxFace f = eb.GetFace(c, q);
          sum += m.GetOutwardSurface(c, q) * feu[f];
        }
        fcg[c] = sum / m.GetVolume(c);
      } else { // Type::cut
        Vect sum = eb.GetSurface(c) * feu[c];
        for (auto q : eb.Nci(c)) {
          const IdxFace f = eb.GetFace(c, q);
          sum += eb.GetOutwardSurface(c, q) * feu[f];
        }
        fcg[c] = sum / eb.GetVolume(c);
      }
    }
    return fcg;
  }
  // Gradient from linear fit to face centers.
  // feu: field on embedded boundaries [a]
  // Returns:
  // gradient on cells [a]
  static FieldCell<Vect> GradientLinearFit(
      const FieldEmbed<Scal>& feu, const EB& eb) {
    auto& m = eb.GetMesh();
    FieldCell<Vect> fcg(eb, Vect(0));
    for (auto c : eb.AllCells()) {
      if (eb.GetType(c) == Type::regular) {
        Vect sum(0);
        for (auto q : m.Nci(c)) {
          const IdxFace f = eb.GetFace(c, q);
          sum += m.GetOutwardSurface(c, q) * feu[f];
        }
        fcg[c] = sum / m.GetVolume(c);
      } else { // Type::cut
        std::vector<Vect> xx;
        std::vector<Scal> uu;
        xx.push_back(eb.GetFaceCenter(c));
        uu.push_back(feu[c]);
        for (auto q : eb.Nci(c)) {
          const IdxFace f = eb.GetFace(c, q);
          xx.push_back(eb.GetFaceCenter(f));
          uu.push_back(feu[f]);
        }
        auto p = ULinear<Scal>::FitLinear(xx, uu);
        fcg[c] = p.first;
      }
    }
    return fcg;
  }
  template <class T>
  static FieldCell<T> AverageCutCells(const FieldCell<T>& fcu, const EB& eb) {
    FieldCell<T> fcr = fcu;
    for (auto c : eb.Cells()) {
      if (eb.GetType(c) == Type::cut) {
        const Scal v = eb.GetVolume(c);
        T sum = fcu[c] * v;
        Scal sumv = v;
        for (IdxCell cn : eb.Stencil(c)) {
          const Scal vn = eb.GetVolume(cn);
          sum += fcu[cn] * vn;
          sumv += vn;
        }
        fcr[c] = sum / sumv;
      }
    }
    return fcr;
  }
  template <class T>
  static FieldCell<T> RedistributeCutCells(
      const FieldCell<T>& fcu, const EB& eb) {
    auto& m = eb.GetMesh();
    FieldCell<T> fcr = fcu;
    for (auto c : eb.Cells()) {
      const Scal v0 = m.GetVolume(c);
      const Scal v = eb.GetVolume(c);
      // excess quantity
      const T du = fcu[c] * (1 - v / v0);
      // subtract from current cell
      fcr[c] -= du;
      // add from neighbor cells proportional to their volume
      for (auto cn : eb.Stencil(c)) {
        if (c != cn) {
          const Scal vn = eb.GetVolume(cn);
          // excess quantity in cell cn
          const T dun = fcu[cn] * (1 - vn / v0);
          fcr[c] += dun * (v / (eb.GetVolumeStencilSum(cn) - vn));
        }
      }
    }
    return fcr;
  }
  // Gradient on faces with limited denominator.
  // fcu: field [a]
  // bc: boundary conditions type, 0: value, 1: grad
  // bcv: value or normal gradient (grad dot GetNormal)
  // Returns:
  // grad dot GetNormal on embedded boundaries [s]
  template <class T>
  static FieldEmbed<T> Gradient(
      const FieldCell<T>& fcu, const MapCondFace& mfc, size_t bc, T bcv,
      const EB& eb) {
    FieldEmbed<T> feg(eb, T(0));
    for (auto f : eb.SuFaces()) {
      const IdxCell cm = eb.GetCell(f, 0);
      const IdxCell cp = eb.GetCell(f, 1);
      const Scal dn = eb.ClipGradDenom(
          (eb.GetNormal(f)).dot(eb.GetCellCenter(cp) - eb.GetCellCenter(cm)));
      feg[f] = (fcu[cp] - fcu[cm]) / dn;
    }
    for (auto c : eb.SuCFaces()) {
      if (bc == 0) {
        const Scal dn = eb.ClipGradDenom(eb.GetFaceOffset(c));
        feg[c] = (bcv - fcu[c]) / dn;
      } else if (bc == 1) {
        feg[c] = bcv;
      } else {
        throw std::runtime_error("Gradient: unknown bc=" + std::to_string(bc));
      }
    }
    GradientB(fcu, mfc, eb.GetMesh(), feg.GetFieldFace());
    return feg;
  }
  // Gradient from iterative interpolation.
  // fcu: field [a]
  // bc: boundary conditions type, 0: value, 1: grad
  // bcv: value or normal gradient (grad dot GetNormal)
  // Returns:
  // grad dot GetNormal on embedded boundaries [s]
  template <class T>
  static FieldEmbed<T> GradientIterative(
      const FieldCell<T>& fcu, const EB& eb) {
    auto feu = InterpolateIterative(fcu, eb);
    auto fcg = GradientLinearFit(feu, eb);
    FieldEmbed<T> feg(eb, T(0));
    for (auto f : eb.SuFaces()) {
      const IdxCell cm = eb.GetCell(f, 0);
      const IdxCell cp = eb.GetCell(f, 1);
      if (eb.GetType(cm) == Type::regular && eb.GetType(cp) == Type::regular) {
        feg[f] = (fcu[cp] - fcu[cm]) / eb.GetCellSize()[0];
      } else {
        feg[f] = (fcg[cm] + fcg[cp]).dot(eb.GetNormal(f)) * 0.5;
      }
    }
    for (auto c : eb.SuCFaces()) {
      feg[c] = fcg[c].dot(eb.GetNormal(c));
    }
    return feg;
  }
  // Updates flux in cut faces using bilinear interpolation from regular faces
  // (Schwartz,2006).
  // feu: field [a]
  // bc: boundary conditions type, 0: value, 1: grad
  // bcv: value or normal gradient (grad dot GetNormal)
  // Returns:
  // grad dot GetNormal on embedded boundaries [s]
  template <class T>
  static FieldFace<T> InterpolateBilinearFaces(
      const FieldFace<T>& ffu, const EB& eb) {
    auto& m = eb.GetMesh();
    FieldFace<T> ffb(eb, T(0));
    for (auto f : eb.SuFaces()) {
      const IdxCell cm = eb.GetCell(f, 0);
      const IdxCell cp = eb.GetCell(f, 1);
      const Scal h = eb.GetCellSize()[0];
      if (eb.GetType(cm) == Type::regular && eb.GetType(cp) == Type::regular) {
        ffb[f] = ffu[f];
      } else {
        const size_t qz = eb.GetNci(cm, f); // f == GetFace(c, qz)
        //                                             //
        //            ---------------------            //
        //            |         |         |            //
        //            |   c01   |   c11   |            //
        //            |         |         |            //
        //            |\--------|---------|            //
        //            | \       |         |            //
        //            |  \ c00  |   c10   |            //
        //            |   \     |         |            //
        //            -----\---------------            //
        //                                             //
        const size_t dz = qz / 2; // face direction
        const size_t dx = (dz + 1) % dim; // other directions
        const size_t dy = (dz + 2) % dim; // other directions
        const Vect n = eb.GetNormal(cm);
        const size_t qx = dx * 2 + (n[dx] < 0 ? 1 : 0);
        const size_t qy = dy * 2 + (n[dy] < 0 ? 1 : 0);
        const IdxCell c00 = cm;
        const IdxCell c10 = eb.GetCell(c00, qx);
        const IdxCell c01 = eb.GetCell(c00, qy);
        const IdxCell c11 = eb.GetCell(c10, qy);
        const IdxFace f00 = eb.GetFace(c00, qz);
        const IdxFace f10 = eb.GetFace(c10, qz);
        const IdxFace f01 = eb.GetFace(c01, qz);
        const IdxFace f11 = eb.GetFace(c11, qz);

        assert(f == f00);
        const Scal tx =
            1 - std::abs(eb.GetFaceCenter(f)[dx] - m.GetCenter(f)[dx]) / h;
        const Scal ty =
            1 - std::abs(eb.GetFaceCenter(f)[dy] - m.GetCenter(f)[dy]) / h;
        const T fy0 = ffu[f00] * tx + ffu[f10] * (1 - tx);
        const T fy1 = ffu[f01] * tx + ffu[f11] * (1 - tx);
        ffb[f] = fy0 * ty + fy1 * (1 - ty);
      }
    }
    return ffb;
  }
  // Interpolates field from cells to faces, bilinear in cut faces.
  // fcu: field on cells [s]
  // Returns:
  // ffu: field on faces [i]
  template <class T>
  static FieldFace<T> InterpolateBilinear(
      const FieldCell<T>& fcu, const EB& eb) {
    FieldFace<T> ffu(eb, T(0));
    for (auto f : eb.SuFaces()) {
      const IdxCell cm = eb.GetCell(f, 0);
      const IdxCell cp = eb.GetCell(f, 1);
      ffu[f] = (fcu[cp] + fcu[cm]) * 0.5;
    }
    return InterpolateBilinearFaces(ffu, eb);
  }
  // Interpolates field from cells to embed faces.
  // fcu: field on cells [s]
  // mcu: value of on embed faces
  // Returns:
  // feu: field on embed faces [i]
  template <class T>
  static void InterpolateEmbedFaces(
      const FieldCell<T>& fcu, size_t bc, const MapCell<T>& mcu,
      FieldEmbed<T>& feu, const EB& eb) {
    auto& m = eb.GetMesh();
    if (bc == 0) {
      for (auto c : eb.SuCFaces()) {
        feu[c] = mcu.at(c);
      }
    } else if (bc == 1) {
      for (auto c : eb.SuCFaces()) {
        std::vector<Vect> xx;
        std::vector<T> uu;
        for (auto cn : eb.Stencil(c)) {
          xx.push_back(m.GetCenter(cn));
          uu.push_back(fcu[cn]);
        }
        auto p = ULinear<Scal>::FitLinear(xx, uu);
        const Scal h = eb.GetCellSize()[0];
        const Vect xc = eb.GetFaceCenter(c) - eb.GetNormal(c) * h;
        const T uc = ULinear<Scal>::EvalLinear(p, xc);
        feu[c] = uc + mcu.at(c) * h;
      }
    } else {
      throw std::runtime_error("not implemented");
    }
  }
  // Interpolates field from cells to embed faces.
  // fcu: field on cells [s]
  // mcu: value of on embed faces
  // Returns:
  // feu: field on embed faces [i]
  template <class T>
  static void InterpolateUpwindEmbedFaces(
      const FieldCell<T>& fcu, size_t bc, const MapCell<Scal>& mcu,
      const FieldEmbed<Scal>& fev, FieldEmbed<T>& feu, const EB& eb) {
    auto& m = eb.GetMesh();
    if (bc == 0) {
      for (auto c : eb.SuCFaces()) {
        if (fev[c] > 0) {
          feu[c] = fcu[c];
        } else if (fev[c] <= 0) {
          feu[c] = mcu.at(c);
        }
      }
    } else if (bc == 1) {
      for (auto c : eb.SuCFaces()) {
        if (fev[c] > 0) {
          feu[c] = fcu[c];
        } else if (fev[c] <= 0) {
          std::vector<Vect> xx;
          std::vector<Scal> uu;
          for (auto cn : eb.Stencil(c)) {
            xx.push_back(m.GetCenter(cn));
            uu.push_back(fcu[cn]);
          }
          auto p = ULinear<Scal>::FitLinear(xx, uu);
          const Scal h = eb.GetCellSize()[0];
          const Vect xc = eb.GetFaceCenter(c) - eb.GetNormal(c) * h;
          const Scal uc = ULinear<Scal>::EvalLinear(p, xc);
          feu[c] = uc + mcu.at(c) * h;
        }
      }
    } else {
      throw std::runtime_error("not implemented");
    }
  }
  // Interpolates field from cells to combined field, bilinear in cut faces.
  // fcu: field on cells [s]
  // mcu: value of on embed faces
  // Returns:
  // feu: field on embed faces [i]
  template <class T>
  static FieldEmbed<T> InterpolateBilinear(
      const FieldCell<T>& fcu, size_t bc, const MapCell<T>& mcu, const EB& eb) {
    FieldEmbed<T> feu(eb, T(0));
    feu.GetFieldFace() = InterpolateBilinear(fcu, eb);
    InterpolateEmbedFaces(fcu, bc, mcu, feu, eb);
    return feu;
  }
  template <class T>
  static FieldEmbed<T> InterpolateBilinear(
      const FieldCell<T>& fcu, size_t bc, T bcv, const EB& eb) {
    MapCell<T> mcu;
    for (auto c : eb.SuCFaces()) {
      mcu[c] = bcv;
    }
    return InterpolateBilinear(fcu, bc, mcu, eb);
  }
  template <class T>
  static FieldEmbed<T> InterpolateBilinear(
      const FieldCell<T>& fcu, const MapCondFace& mfc, size_t bc, T bcv,
      const EB& eb) {
    MapCell<T> mcu;
    for (auto c : eb.SuCFaces()) {
      mcu[c] = bcv;
    }
    auto feu = InterpolateBilinear(fcu, bc, mcu, eb);
    InterpolateB(fcu, mfc, feu.GetFieldFace(), eb.GetMesh());
    return feu;
  }
  // Gradient with bilinear interpolation in cut faces
  // and linear fit in embed faces.
  // fcu: field [a]
  // Returns:
  // grad dot GetNormal on faces [s]
  template <class T>
  static FieldFace<T> GradientBilinear(const FieldCell<T>& fcu, const EB& eb) {
    FieldFace<T> ffg(eb, T(0));
    const Scal h = eb.GetCellSize()[0];
    for (auto f : eb.SuFaces()) {
      const IdxCell cm = eb.GetCell(f, 0);
      const IdxCell cp = eb.GetCell(f, 1);
      ffg[f] = (fcu[cp] - fcu[cm]) / h;
    }
    return InterpolateBilinearFaces(ffg, eb);
  }
  // Gradient with bilinear interpolation in cut faces
  // and linear fit in embed faces.
  // fcu: field [a]
  // mcu: value (bc=0) or derivative (bc=1)
  // Returns:
  // grad dot GetNormal on embedded boundaries [s]
  template <class T>
  static FieldEmbed<T> GradientBilinear(
      const FieldCell<T>& fcu, size_t bc, const MapCell<T>& mcu, const EB& eb) {
    auto& m = eb.GetMesh();
    FieldEmbed<T> feg(eb, T(0));
    feg.GetFieldFace() = GradientBilinear(fcu, eb);
    if (bc == 0) {
      for (auto c : eb.SuCFaces()) {
        std::vector<Vect> xx;
        std::vector<T> uu;
        for (auto cn : eb.Stencil(c)) {
          xx.push_back(m.GetCenter(cn));
          uu.push_back(fcu[cn]);
        }
        auto p = ULinear<Scal>::FitLinear(xx, uu);
        const T ub = mcu.at(c);
        const Scal h = eb.GetCellSize()[0];
        const Vect xc = eb.GetFaceCenter(c) - eb.GetNormal(c) * h;
        const T uc = ULinear<Scal>::EvalLinear(p, xc);
        feg[c] = (ub - uc) / h;
      }
    } else if (bc == 1) {
      for (auto c : eb.SuCFaces()) {
        feg[c] = mcu.at(c);
      }
    } else {
      throw std::runtime_error("not implemented");
    }
    return feg;
  }
  template <class T>
  static FieldEmbed<T> GradientBilinear(
      const FieldCell<T>& fcu, const MapCondFace& mfc, size_t bc, T bcv,
      const EB& eb) {
    MapCell<T> mcu;
    for (auto c : eb.SuCFaces()) {
      mcu[c] = bcv;
    }
    auto feg = GradientBilinear(fcu, bc, mcu, eb);
    GradientB(fcu, mfc, eb.GetMesh(), feg.GetFieldFace());
    return feg;
  }
};
