// Created by Petr Karnakov on 03.02.2020
// Copyright 2020 ETH Zurich

#undef NDEBUG
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>

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
  virtual ~Func() = default;
  virtual std::function<T(Vect)> operator()() const = 0;
  virtual std::function<T(Vect)> Dx(size_t d) const = 0;
  virtual std::function<T(Vect)> Dxx(size_t d, size_t dd) const = 0;
};

template <class T>
class Quad : public Func<T> {
 public:
  Quad() = default;
  virtual ~Quad() = default;
  std::function<T(Vect)> operator()() const override {
    return [](Vect v) {
      const T x(v[0]), y(v[1]), z(v[2]);
      return x * x + y * y + z * z;
    };
  }
  std::function<T(Vect)> Dx(size_t d) const override {
    return [d](Vect v) { return T(2 * v[d]); };
  }
  std::function<T(Vect)> Dxx(size_t d, size_t dd) const override {
    return [d, dd](Vect) { return T(2 * (d == dd)); };
  }
};

class Sin : public Func<Scal> {
 public:
  Sin(size_t seed = 0) : om_(std::sin(seed + 1) * 10) {}
  virtual ~Sin() = default;
  std::function<Scal(Vect)> operator()() const override {
    return [this](Vect v) {
      v *= om_;
      const Scal x(v[0]), y(v[1]), z(v[2]);
      return std::sin(x) * std::sin(y) * std::sin(z);
    };
  }
  std::function<Scal(Vect)> Dx(size_t d) const override {
    return [d, this](Vect v) {
      v *= om_;
      Scal r = 1;
      for (auto q : Dirs) {
        r *= (q == d ? om_[q] * std::cos(v[q]) : std::sin(v[q]));
      }
      return r;
    };
  }
  std::function<Scal(Vect)> Dxx(size_t d, size_t dd) const override {
    return [d, dd, this](Vect v) {
      v *= om_;
      Scal r = 1;
      if (d == dd) {
        for (auto q : Dirs) {
          r *= (q == d ? -sqr(om_[q]) * std::sin(v[q]) : std::sin(v[q]));
        }
      } else {
        for (auto q : Dirs) {
          r *= (q == d || q == dd ? om_[q] * std::cos(v[q]) : std::sin(v[q]));
        }
      }
      return r;
    };
  }

 private:
  const Vect om_;
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
  const MIdx size(32);
  Rect<Vect> dom(Vect(0), Vect(h) * Vect(size));
  return InitUniformMesh<M>(dom, MIdx(0), size, 2, true, true, size, 0);
}

template <class Idx, class M>
Scal GetNorm1(const GField<Scal, Idx>& u, const M& m) {
  Scal sum = 0;
  Scal sumv = 0;
  for (auto i : m.template GetAll<Idx>()) {
    sum += std::abs(u[i]) * m.GetVolume(i);
    sumv += m.GetVolume(i);
  }
  return sum / sumv;
}

template <class Idx>
Scal GetNorm2(const GField<Scal, Idx>& u, const M& m) {
  Scal sum = 0;
  Scal sumv = 0;
  for (auto i : m.template GetAll<Idx>()) {
    sum += sqr(u[i]) * m.GetVolume(i);
    sumv += m.GetVolume(i);
  }
  return std::sqrt(sum / sumv);
}

template <class T, class Idx>
GField<T, Idx> Subtract(
    const GField<T, Idx>& u0, const GField<T, Idx>& u1, const M& m) {
  GField<Scal, Idx> d(m);
  for (auto i : m.template GetAll<Idx>()) {
    d[i] = u0[i] - u1[i];
  }
  return d;
}

int main() {
  const Scal h0 = 1. / 32;
  for (auto h : {h0, h0 / 2, h0 / 4}) {
    auto m = GetMesh(h);
    // Quad<Scal> u;
    Sin u;
    auto fcu = GetField<Scal, IdxCell>(u(), m);
    auto ffu = GetField<Scal, IdxFace>(u(), m);
    auto fcui = Average(ffu, m);
    std::cout << GetNorm2(Subtract(fcu, fcui, m), m) << std::endl;
  }
}
