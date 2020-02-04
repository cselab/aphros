// Created by Petr Karnakov on 03.02.2020
// Copyright 2020 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <random>
#include <sstream>
#include <typeinfo>

#include "geom/mesh.h"
#include "solver/approx.h"
#include "solver/cond.h"
#include "solver/solver.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;
auto Dirs = GRange<size_t>(3);

template <class T>
class Func {
 public:
  explicit Func(int seed) : gen_(seed) {}
  virtual ~Func() = default;
  virtual std::function<T(Vect)> operator()() const = 0;
  virtual std::function<T(Vect)> Dx(size_t d) const = 0;
  virtual std::function<T(Vect)> Dxx(size_t d, size_t dd) const = 0;

 protected:
  Scal G(Scal min = 0, Scal max = 1) {
    return std::uniform_real_distribution<Scal>(min, max)(gen_);
  }
  std::default_random_engine gen_;
};

class Quad : public Func<Scal> {
 public:
  explicit Quad(int seed) : Func<Scal>(seed) {}
  virtual ~Quad() = default;
  std::function<Scal(Vect)> operator()() const override {
    return [](Vect v) {
      v += Vect(1);
      const Scal x(v[0]), y(v[1]), z(v[2]);
      return x * x + y * y + z * z;
    };
  }
  std::function<Scal(Vect)> Dx(size_t d) const override {
    return [d](Vect v) {
      v += Vect(1);
      return 2 * v[d];
    };
  }
  std::function<Scal(Vect)> Dxx(size_t d, size_t dd) const override {
    return [d, dd](Vect) { return 2 * (d == dd); };
  }
};

class Sin : public Func<Scal> {
 public:
  explicit Sin(int seed)
      : Func<Scal>(seed)
      , freq_(Vect(G(), G(), G()) + Vect(10))
      , phase_(Vect(0.5) + Vect(G(), G(), G()) * 10) {}
  virtual ~Sin() = default;
  std::function<Scal(Vect)> operator()() const override {
    return [this](Vect v) {
      v = freq_ * v + phase_;
      const Scal x(v[0]), y(v[1]), z(v[2]);
      return std::sin(x) * std::sin(y) * std::sin(z);
    };
  }
  std::function<Scal(Vect)> Dx(size_t d) const override {
    return [d, this](Vect v) {
      v = freq_ * v + phase_;
      Scal r = 1;
      for (auto q : Dirs) {
        r *= (q == d ? freq_[q] * std::cos(v[q]) : std::sin(v[q]));
      }
      return r;
    };
  }
  std::function<Scal(Vect)> Dxx(size_t d, size_t dd) const override {
    return [d, dd, this](Vect v) {
      v = freq_ * v + phase_;
      Scal r = 1;
      if (d == dd) {
        for (auto q : Dirs) {
          r *= (q == d ? -sqr(freq_[q]) * std::sin(v[q]) : std::sin(v[q]));
        }
      } else {
        for (auto q : Dirs) {
          r *= (q == d || q == dd ? freq_[q] * std::cos(v[q]) : std::sin(v[q]));
        }
      }
      return r;
    };
  }

 private:
  const Vect freq_;
  const Vect phase_;
};

std::function<Vect(Vect)> Gradient(const Func<Scal>& u) {
  return [&u](Vect x) { return Vect(u.Dx(0)(x), u.Dx(1)(x), u.Dx(2)(x)); };
}

template <class T, class Idx>
GField<T, Idx> GetField(const std::function<T(Vect)>& u, const M& m) {
  GField<T, Idx> r(m);
  for (auto i : m.template GetAll<Idx>()) {
    r[i] = u(m.GetCenter(i));
  }
  return r;
}

// h: cell size
M GetMesh(Scal h) {
  const MIdx size(16);
  Rect<Vect> dom(Vect(0), Vect(h) * Vect(size));
  return InitUniformMesh<M>(dom, MIdx(0), size, 2, true, true, size, 0);
}

template <class Idx, class M>
Scal Norm1(const GField<Scal, Idx>& u, const M& m) {
  Scal sum = 0;
  Scal sumv = 0;
  for (auto i : m.template GetIn<Idx>()) {
    sum += std::abs(u[i]);
    sumv += 1;
  }
  return sum / sumv;
}

template <class Idx>
Scal Norm2(const GField<Scal, Idx>& u, const M& m) {
  Scal sum = 0;
  Scal sumv = 0;
  for (auto i : m.template GetIn<Idx>()) {
    sum += sqr(u[i]);
    sumv += 1;
  }
  return std::sqrt(sum / sumv);
}

template <class T, class Idx>
GField<T, Idx> Add(
    const GField<T, Idx>& u0, const GField<T, Idx>& u1, const M& m) {
  GField<Scal, Idx> d(m);
  for (auto i : m.template GetAll<Idx>()) {
    d[i] = u0[i] + u1[i];
  }
  return d;
}

template <class T, class Idx>
GField<T, Idx> Subtract(
    const GField<T, Idx>& u0, const GField<T, Idx>& u1, const M& m) {
  GField<Scal, Idx> r(m);
  for (auto i : m.template GetAll<Idx>()) {
    r[i] = u0[i] - u1[i];
  }
  return r;
}

template <class T, class Idx>
GField<T, Idx> Mul(const GField<T, Idx>& u, Scal k, const M& m) {
  GField<Scal, Idx> r(m);
  for (auto i : m.template GetAll<Idx>()) {
    r[i] = u[i] * k;
  }
  return r;
}

template <class T, class Idx>
GField<T, Idx> Abs(const GField<T, Idx>& u, const M& m) {
  GField<Scal, Idx> r(m);
  for (auto i : m.template GetAll<Idx>()) {
    r[i] = std::abs(u[i]);
  }
  return r;
}

// Evaluates error on meshes and estimates convergence order.
// F: function to evaluate, derived from Func, constructable from int
// func: function to evaluate
// estimator: estimator of field
// exact: exact field
// Returns:
// tuple(error0, error1, order)
template <class F, class T, class Idx>
std::tuple<Scal, Scal, Scal> EstimateOrder(
    std::function<GField<Scal, Idx>(const Func<T>&, const M&)> estimator,
    std::function<GField<Scal, Idx>(const Func<T>&, const M&)> exact) {
  using Field = GField<T, Idx>;
  const Scal h0 = 1e-3;
  const size_t nsamp = 10;
  const std::vector<Scal> hh = {h0, h0 * 0.5};
  std::vector<Scal> ee(hh.size());
  for (size_t i = 0; i < hh.size(); ++i) {
    const auto m = GetMesh(hh[i]);
    Field avg(m, 0);
    for (size_t seed = 0; seed < nsamp; ++seed) {
      const F func(seed);
      auto f0 = estimator(func, m);
      auto f1 = exact(func, m);
      avg = Add(avg, Abs(Subtract(f0, f1, m), m), m);
    }
    avg = Mul(avg, 1. / nsamp, m);
    ee[i] = Norm1(avg, m);
  }
  return std::make_tuple(
      ee[0], ee[1], std::log(ee[0] / ee[1]) / std::log(hh[0] / hh[1]));
}

template <class F, class T, class Idx>
void PrintEstimateOrder(
    std::function<GField<T, Idx>(const Func<Scal>&, const M&)> estimator,
    std::function<GField<T, Idx>(const Func<Scal>&, const M&)> exact) {
  std::cout << "func=" << typeid(F).name() << std::endl;
  auto t = EstimateOrder<F, T, Idx>(estimator, exact);
  printf(
      "order=%5.3f   coarse=%.5e   fine=%.5e\n", std::get<2>(t), std::get<0>(t),
      std::get<1>(t));
}

template <class T, class Idx>
void VaryFunc(
    std::function<GField<T, Idx>(const Func<Scal>&, const M&)> estimator,
    std::function<GField<T, Idx>(const Func<Scal>&, const M&)> exact) {
  PrintEstimateOrder<Quad, Scal, Idx>(estimator, exact);
  PrintEstimateOrder<Sin, Scal, Idx>(estimator, exact);
}

int main() {
  std::cout << "\nGradientI() returning FieldFace<Scal>" << std::endl;
  VaryFunc<Scal, IdxFace>(
      [](const Func<Scal>& func, const M& m) {
        const auto u = GetField<Scal, IdxCell>(func(), m);
        FieldFace<Scal> g(m);
        GradientI(u, m, g);
        return g;
      },
      [](const Func<Scal>& func, const M& m) {
        FieldFace<Scal> g(m, 0);
        for (auto f : m.AllFaces()) {
          const auto x = m.GetCenter(f);
          g[f] = Gradient(func)(x).dot(m.GetNormal(f));
        }
        return g;
      });

  std::cout << "\nAverage() returning FieldCell<Scal> " << std::endl;
  VaryFunc<Scal, IdxCell>(
      [](const Func<Scal>& func, const M& m) {
        const auto ffu = GetField<Scal, IdxFace>(func(), m);
        const auto fcu = Average(ffu, m);
        return fcu;
      },
      [](const Func<Scal>& func, const M& m) {
        return GetField<Scal, IdxCell>(func(), m);
      });
}
