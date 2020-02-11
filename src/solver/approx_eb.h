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

#if 0

template <class M>
struct UEmbed {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  constexpr size_t dim = 3;

  template <class Scal, size_t N>
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
  template <class Scal, size_t N>
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
  std::pair<Vect, Scal> FitLinear(
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
  template <class Scal>
  std::pair<generic::Vect<Vect, dim>, Vect> FitLinear(
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

  template <class T, class Scal>
  T EvalLinear(
      const std::pair<generic::Vect<T, dim>, T>& p,
      const generic::Vect<Scal, dim>& x) {
    auto& g = p.first;
    auto& u0 = p.second;
    return g[0] * x[0] + g[1] * x[1] + g[2] * x[2] + u0;
  }

  template <class Scal>
  struct CondEmbed {
    enum class Type { value, gradient };
    Scal u; // value or normal gradient
  };

  FieldNode<Scal> InitEmbed(const M& m, const Vars& var, bool verb) {
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
};

#endif
