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

#include "dump/dump.h"
#include "geom/mesh.h"
#include "solver/approx.h"
#include "solver/approx_eb.h"
#include "solver/cond.h"
#include "solver/embed.h"
#include "solver/solver.h"
#include "util/fluid.h"

using M = MeshCartesian<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;
auto Dirs = GRange<size_t>(3);
using EB = Embed<M>;
using UEB = UEmbed<M>;

template <class T>
class Func {
 public:
  using Result = T;
  explicit Func(int seed) : gen_(seed) {}
  virtual ~Func() = default;
  virtual std::function<T(Vect)> operator()() const = 0;
  virtual std::function<T(Vect)> Dx(size_t d) const = 0;
  virtual std::function<T(Vect)> Dxx(size_t d, size_t dd) const = 0;
  static std::string GetName() {
    return "";
  }

 protected:
  Scal G(Scal min = 0, Scal max = 1) {
    return std::uniform_real_distribution<Scal>(min, max)(gen_);
  }
  std::mt19937 gen_;
};

class Quadratic final : public Func<Scal> {
 public:
  explicit Quadratic(int seed)
      : Func<Scal>(seed), origin_(Vect(1) + Vect(G(), G(), G()) * 0.1) {}
  std::function<Scal(Vect)> operator()() const override {
    return [this](Vect v) {
      v -= origin_;
      const Scal x(v[0]), y(v[1]), z(v[2]);
      return x * x + y * y + z * z;
    };
  }
  std::function<Scal(Vect)> Dx(size_t d) const override {
    return [d, this](Vect v) {
      v -= origin_;
      return 2 * v[d];
    };
  }
  std::function<Scal(Vect)> Dxx(size_t d, size_t dd) const override {
    return [d, dd](Vect) { return 2 * (d == dd); };
  }
  static std::string GetName() {
    return "Quadratic";
  }

 private:
  const Vect origin_;
};

class QuadraticBilinear final : public Func<Scal> {
 public:
  const size_t dim = Vect::dim;
  explicit QuadraticBilinear(int seed)
      : Func<Scal>(seed), origin_(Vect(1) + Vect(G(), G(), G()) * 0.1) {}
  std::function<Scal(Vect)> operator()() const override {
    return [this](Vect v) {
      v -= origin_;
      const Scal x(v[0]), y(v[1]), z(v[2]);
      return x * x * y * z + y * y * x * z + z * z * x * y;
    };
  }
  std::function<Scal(Vect)> Dx(size_t d) const override {
    return [d, this](Vect v) {
      v -= origin_;
      const size_t dx = d;
      const size_t dy = (d + 1) % dim;
      const size_t dz = (d + 2) % dim;
      const Scal x(v[dx]), y(v[dy]), z(v[dz]);
      return 2 * x * y * z + y * y * z + z * z * y;
    };
  }
  std::function<Scal(Vect)> Dxx(size_t d, size_t dd) const override {
    return [d, dd, this](Vect v) {
      v -= origin_;
      const size_t dx = d;
      size_t dy = (d + 1) % dim;
      size_t dz = (d + 2) % dim;
      if (d == dd) {
        const Scal y(v[dy]), z(v[dz]);
        return 2 * y * z;
      }
      if (dy != dd) {
        std::swap(dy, dz);
      }
      const Scal x(v[dx]), y(v[dy]), z(v[dz]);
      return 2 * x * z + 2 * y * z + z * z;
    };
  }
  static std::string GetName() {
    return "QuadraticBilinear";
  }

 private:
  const Vect origin_;
};

class Trilinear final : public Func<Scal> {
 public:
  const size_t dim = Vect::dim;
  explicit Trilinear(int seed)
      : Func<Scal>(seed), origin_(Vect(1) + Vect(G(), G(), G()) * 0.1) {}
  std::function<Scal(Vect)> operator()() const override {
    return [this](Vect v) {
      v -= origin_;
      const Scal x(v[0]), y(v[1]), z(v[2]);
      return x * y * z;
    };
  }
  std::function<Scal(Vect)> Dx(size_t d) const override {
    return [d, this](Vect v) {
      v -= origin_;
      const size_t dy = (d + 1) % dim;
      const size_t dz = (d + 2) % dim;
      const Scal y(v[dy]), z(v[dz]);
      return y * z;
    };
  }
  std::function<Scal(Vect)> Dxx(size_t d, size_t dd) const override {
    return [d, dd, this](Vect v) {
      v -= origin_;
      size_t dy = (d + 1) % dim;
      size_t dz = (d + 2) % dim;
      if (d == dd) {
        return 0.;
      }
      if (dy != dd) {
        std::swap(dy, dz);
      }
      return v[dz];
    };
  }
  static std::string GetName() {
    return "Trilinear";
  }

 private:
  const Vect origin_;
};

class Linear final : public Func<Scal> {
 public:
  explicit Linear(int seed)
      : Func<Scal>(seed), origin_(Vect(1) + Vect(G(), G(), G()) * 0.1) {}
  std::function<Scal(Vect)> operator()() const override {
    return [this](Vect v) {
      v -= origin_;
      return v.dot(Vect(1));
    };
  }
  std::function<Scal(Vect)> Dx(size_t) const override {
    return [this](Vect) { return 1; };
  }
  std::function<Scal(Vect)> Dxx(size_t, size_t) const override {
    return [](Vect) { return 0; };
  }
  static std::string GetName() {
    return "Linear";
  }

 private:
  const Vect origin_;
};

class Const final : public Func<Scal> {
 public:
  explicit Const(int seed) : Func<Scal>(seed), const_(1 + G() * 0.1) {}
  std::function<Scal(Vect)> operator()() const override {
    return [this](Vect) { return const_; };
  }
  std::function<Scal(Vect)> Dx(size_t) const override {
    return [this](Vect) { return 0; };
  }
  std::function<Scal(Vect)> Dxx(size_t, size_t) const override {
    return [](Vect) { return 0; };
  }
  static std::string GetName() {
    return "Const";
  }

 private:
  const Scal const_;
};

class Sine final : public Func<Scal> {
 public:
  explicit Sine(int seed)
      : Func<Scal>(seed)
      , freq_(Vect(G(), G(), G()) + Vect(10))
      , phase_(Vect(0.5) + Vect(G(), G(), G()) * 10) {}
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
  static std::string GetName() {
    return "Sine";
  }

 private:
  const Vect freq_;
  const Vect phase_;
};

std::function<Vect(Vect)> Gradient(const Func<Scal>& u) {
  return [&u](Vect x) { return Vect(u.Dx(0)(x), u.Dx(1)(x), u.Dx(2)(x)); };
}

template <class F>
class FuncVect final : public Func<Vect> {
 public:
  FuncVect(int seed) : Func<Vect>(seed), dr_(Vect::dim) {
    for (size_t d = 0; d < Vect::dim; ++d) {
      vfunc_.emplace_back(int(G() * 100));
    }
  };
  std::function<Vect(Vect)> operator()() const override {
    return [this](Vect v) {
      Vect r;
      for (auto i : dr_) {
        r[i] = vfunc_[i]()(v);
      }
      return r;
    };
  }
  std::function<Vect(Vect)> Dx(size_t d) const override {
    return [d, this](Vect v) {
      Vect r;
      for (auto i : dr_) {
        r[i] = vfunc_[i].Dx(d)(v);
      }
      return r;
    };
  }
  std::function<Vect(Vect)> Dxx(size_t d, size_t dd) const override {
    return [d, dd, this](Vect v) {
      Vect r;
      for (auto i : dr_) {
        r[i] = vfunc_[i].Dxx(d, dd)(v);
      }
      return r;
    };
  }
  static std::string GetName() {
    return "Vect" + F::GetName();
  }

 private:
  std::vector<F> vfunc_;
  const GRange<size_t> dr_;
};

namespace eval_impl_1 {
template <class Field, class MEB>
struct EvalHelper {
  using Value = typename Field::Value;
  Field operator()(const std::function<Value(Vect)>& u, const MEB& meb);
};

template <class Field>
struct EvalHelper<Field, M> {
  using Value = typename Field::Value;
  using Idx = typename Field::Idx;
  Field operator()(const std::function<Value(Vect)>& u, const M& m) {
    Field r(m, Value(0));
    for (auto i : m.template GetRangeAll<Idx>()) {
      r[i] = u(m.GetCenter(i));
    }
    return r;
  }
};

template <class T>
struct EvalHelper<FieldCell<T>, EB> {
  using Value = T;
  using Field = FieldCell<Value>;
  Field operator()(const std::function<Value(Vect)>& u, const EB& eb) {
    Field r(eb, Value(0));
    for (auto c : eb.Cells()) {
      r[c] = u(eb.GetCellCenter(c));
    }
    return r;
  }
};

template <class T>
struct EvalHelper<FieldEmbed<T>, EB> {
  using Value = T;
  using Field = FieldEmbed<Value>;
  Field operator()(const std::function<Value(Vect)>& u, const EB& eb) {
    Field r(eb, Value(0));
    for (auto c : eb.CFaces()) {
      r[c] = u(eb.GetFaceCenter(c));
    }
    for (auto f : eb.Faces()) {
      r[f] = u(eb.GetFaceCenter(f));
    }
    return r;
  }
};
} // namespace eval_impl_1

template <class Field, class MEB>
Field Eval(
    const std::function<typename Field::Value(Vect)>& u, const MEB& meb) {
  return eval_impl_1::EvalHelper<Field, MEB>()(u, meb);
}

namespace eval_impl_2 {

template <class Field, class MEB>
struct EvalHelper {
  using Value = typename Field::Value;
  Field operator()(
      const std::function<Value(Value)>& func, const Field& u0, const MEB& meb);
};
template <class Field>
struct EvalHelper<Field, M> {
  using Value = typename Field::Value;
  using Idx = typename Field::Idx;
  Field operator()(
      const std::function<Value(Value)>& func, const Field& u0, const M& m) {
    Field r(m, Value(0));
    for (auto i : m.template GetRangeAll<Idx>()) {
      r[i] = func(u0[i]);
    }
    return r;
  }
};
template <class T>
struct EvalHelper<FieldCell<T>, EB> {
  using Value = T;
  using Field = FieldCell<Value>;
  Field operator()(
      const std::function<Value(Value)>& func, const Field& u0, const EB& eb) {
    Field r(eb, Value(0));
    for (auto c : eb.Cells()) {
      r[c] = func(u0[c]);
    }
    return r;
  }
};
template <class T>
struct EvalHelper<FieldEmbed<T>, EB> {
  using Value = T;
  using Field = FieldEmbed<Value>;
  Field operator()(
      const std::function<Value(Value)>& func, const Field& u0, const EB& eb) {
    Field r(eb, Value(0));
    for (auto c : eb.CFaces()) {
      r[c] = func(u0[c]);
    }
    for (auto f : eb.Faces()) {
      r[f] = func(u0[f]);
    }
    return r;
  }
};

} // namespace eval_impl_2

template <class Field, class MEB>
Field Eval(
    const std::function<typename Field::Value(typename Field::Value)>& func,
    const Field& u0, const MEB& meb) {
  return eval_impl_2::EvalHelper<Field, MEB>()(func, u0, meb);
}

namespace eval_impl_3 {

template <class Field, class MEB>
struct EvalHelper {
  using Value = typename Field::Value;
  Field operator()(
      const std::function<Value(Value, Value)>& func, const Field& u0,
      const Field& u1, const MEB& meb);
};
template <class Field>
struct EvalHelper<Field, M> {
  using Value = typename Field::Value;
  using Idx = typename Field::Idx;
  Field operator()(
      const std::function<Value(Value, Value)>& func, const Field& u0,
      const Field& u1, const M& m) {
    Field r(m, Value(0));
    for (auto i : m.template GetRangeAll<Idx>()) {
      r[i] = func(u0[i], u1[i]);
    }
    return r;
  }
};
template <class T>
struct EvalHelper<FieldCell<T>, EB> {
  using Value = T;
  using Field = FieldCell<Value>;
  Field operator()(
      const std::function<Value(Value, Value)>& func, const Field& u0,
      const Field& u1, const EB& eb) {
    Field r(eb, Value(0));
    for (auto c : eb.Cells()) {
      r[c] = func(u0[c], u1[c]);
    }
    return r;
  }
};
template <class T>
struct EvalHelper<FieldEmbed<T>, EB> {
  using Value = T;
  using Field = FieldEmbed<Value>;
  Field operator()(
      const std::function<Value(Value, Value)>& func, const Field& u0,
      const Field& u1, const EB& eb) {
    Field r(eb, Value(0));
    for (auto c : eb.CFaces()) {
      r[c] = func(u0[c], u1[c]);
    }
    for (auto f : eb.Faces()) {
      r[f] = func(u0[f], u1[f]);
    }
    return r;
  }
};

} // namespace eval_impl_3

template <class Field, class MEB>
Field Eval(
    const std::function<typename Field::Value(
        typename Field::Value, typename Field::Value)>& func,
    const Field& u0, const Field& u1, const MEB& meb) {
  return eval_impl_3::EvalHelper<Field, MEB>()(func, u0, u1, meb);
}

// h: cell size
std::unique_ptr<M> CreateMesh(Scal h) {
  const MIdx size(16);
  Rect<Vect> dom(Vect(0), Vect(h) * Vect(size));
  return std::make_unique<M>(MIdx(0), size, dom, 2, true, true, size, 0);
}

// TODO: rename SuCells to Cells(1)
// TODO: rename AllCells to Cells(2)

std::unique_ptr<EB> CreateEmbed(M& m) {
  auto peb = std::make_unique<EB>(m, 0);
  FieldNode<Scal> fnl(m);
  auto block = m.GetInBlockCells().GetSize();
  auto h = m.GetCellSize();
  for (auto n : m.AllNodes()) {
    const auto x = m.GetNode(n) / h / Vect(block);
    auto dx = Vect(0.5) - x;
    // dx[0] = 0;
    fnl[n] = 0.41 - dx.norm();
  }
  do {
    peb->Init(fnl);
  } while (m.Pending());
  return peb;
}

template <class MEB>
std::unique_ptr<MEB> ConvertMesh(std::unique_ptr<M>& peb);

template <>
std::unique_ptr<M> ConvertMesh<M>(std::unique_ptr<M>& pm) {
  return std::unique_ptr<M>(pm.release());
}

template <>
std::unique_ptr<EB> ConvertMesh<EB>(std::unique_ptr<M>& pm) {
  return CreateEmbed(*pm);
}

template <class Field, class MEB>
Field Add(const Field& u0, const Field& u1, const MEB& meb) {
  using Value = typename Field::Value;
  return Eval([](Value a0, Value a1) { return a0 + a1; }, u0, u1, meb);
}

template <class Field, class MEB>
Field Sub(const Field& u0, const Field& u1, const MEB& meb) {
  using Value = typename Field::Value;
  return Eval([](Value a0, Value a1) { return a0 - a1; }, u0, u1, meb);
}

template <class Field, class MEB>
Field Div(const Field& u0, const Field& u1, const MEB& meb) {
  using Value = typename Field::Value;
  return Eval([](Value a0, Value a1) { return a0 / a1; }, u0, u1, meb);
}

template <class Field, class MEB>
Field Mul(const Field& u0, Scal k, const MEB& meb) {
  using Value = typename Field::Value;
  return Eval([k](Value a0) { return a0 * k; }, u0, meb);
}

Scal Abs(Scal a) {
  return std::abs(a);
}

Vect Abs(Vect a) {
  return a.abs();
}

template <class Field, class MEB>
Field Abs(const Field& u0, const MEB& meb) {
  using Value = typename Field::Value;
  return Eval([](Value a0) { return Abs(a0); }, u0, meb);
}

Scal Norm1(Scal a) {
  return std::abs(a);
}

Scal Norm1(Vect a) {
  return a.norm1();
}

template <class T, class Idx>
Scal Norm1(const GField<T, Idx>& u, const M& m) {
  Scal sum = 0;
  Scal sumv = 0;
  for (auto i : m.GetRangeIn<Idx>()) {
    sum += Norm1(u[i]);
    sumv += 1;
  }
  return sum / sumv;
}

template <class T>
Scal Norm1(const FieldCell<T>& u, const EB& eb) {
  Scal sum = 0;
  Scal sumv = 0;
  for (auto i : eb.Cells()) {
    sum += Norm1(u[i]);
    sumv += 1;
  }
  return sum / sumv;
}
template <class T>
Scal Norm1(const FieldEmbed<T>& u, const EB& eb) {
  Scal sum = 0;
  Scal sumv = 0;
  for (auto c : eb.CFaces()) {
    sum += Norm1(u[c]);
    sumv += 1;
  }
  for (auto f : eb.Faces()) {
    sum += Norm1(u[f]);
    sumv += 1;
  }
  return sum / sumv;
}

template <class F, class Field, class MEB>
Field Eval(
    std::function<Field(const Func<typename F::Result>&, const MEB&)> estimator,
    Scal h) {
  auto pm = CreateMesh(h);
  auto pmeb = ConvertMesh<MEB>(pm);
  auto& meb = *pmeb;
  return estimator(F(0), meb);
}

// Computes mean error field of estimator compared to exact values.
// F: function to evaluate, derived from Func, constructable from int
// estimator: estimator of field
// exact: exact field
// h: mesh step
// Returns:
// mean error field
template <class F, class Field, class MEB>
Field GetErrorField(
    std::function<Field(const Func<typename F::Result>&, const MEB&)> estimator,
    std::function<Field(const Func<typename F::Result>&, const MEB&)> exact,
    Scal h) {
  const size_t nsamp = 10;
  auto pm = CreateMesh(h);
  auto pmeb = ConvertMesh<MEB>(pm);
  auto& meb = *pmeb;
  using Value = typename Field::Value;
  Field avg(meb, Value(0));
  for (size_t seed = 0; seed < nsamp; ++seed) {
    const F func(seed);
    auto f0 = estimator(func, meb);
    auto f1 = exact(func, meb);
    avg = Add(avg, Abs(Sub(f0, f1, meb), meb), meb);
  }
  avg = Mul(avg, 1. / nsamp, meb);
  return avg;
}

template <class T>
T Eval(T a, const std::function<Scal(Scal)> func);

template <>
Scal Eval(Scal a, const std::function<Scal(Scal)> func) {
  return func(a);
}

template <>
Vect Eval(Vect v, const std::function<Scal(Scal)> func) {
  for (size_t i = 0; i < v.size(); ++i) {
    v[i] = func(v[i]);
  }
  return v;
}

// Computes field with estimated order of accuracy.
// F: function to evaluate, derived from Func, constructable from int
// estimator: estimator of field
// exact: exact field
// Returns:
// tuple(error0, error1, order)
template <class F, class Field, class MEB>
Field GetOrderField(
    std::function<Field(const Func<typename F::Result>&, const MEB&)> estimator,
    std::function<Field(const Func<typename F::Result>&, const MEB&)> exact) {
  const Scal h0 = 1e-3;
  const Scal h1 = h0 * 0.5;
  auto er0 = GetErrorField<F, Field, MEB>(estimator, exact, h0);
  auto er1 = GetErrorField<F, Field, MEB>(estimator, exact, h1);
  auto pm = CreateMesh(h0);
  auto pmeb = ConvertMesh<MEB>(pm);
  auto& meb = *pmeb;
  using Value = typename Field::Value;
  return Eval(
      [h0, h1](Value e0, Value e1) {
        if (Vect(e1) == Vect(0)) {
          return Value(0);
        }
        return Eval(e0 / e1, [](Scal a) { return std::log(a); }) /
               std::log(h0 / h1);
      },
      er0, er1, meb);
}

/// Computes mean error field of estimator compared to exact values.
/// @param F function to evaluate, derived from Func, constructable from int
/// @param func function to evaluate
/// @param estimator estimator of field
/// @param exact exact field
/// @param h mesh step
/// @return tuple(error0, error1, order)

// Evaluates error on two meshes and estimates the order of accuracy.
// F: function to evaluate, derived from Func, constructable from int
// func: function to evaluate
// estimator: estimator of field
// exact: exact field
// Returns:
// tuple(error0, error1, order)
template <class F, class Field, class MEB>
std::tuple<Scal, Scal, Scal> CalcOrder(
    std::function<Field(const Func<typename F::Result>&, const MEB&)> estimator,
    std::function<Field(const Func<typename F::Result>&, const MEB&)> exact) {
  const Scal h0 = 1e-3;
  const std::vector<Scal> hh = {h0, h0 * 0.5};
  std::vector<Scal> ee(hh.size());
  for (size_t i = 0; i < hh.size(); ++i) {
    auto pm = CreateMesh(hh[i]);
    auto pmeb = ConvertMesh<MEB>(pm);
    ee[i] = Norm1(GetErrorField<F, Field, MEB>(estimator, exact, hh[i]), *pmeb);
  }
  return std::make_tuple(
      ee[0], ee[1], std::log(ee[0] / ee[1]) / std::log(hh[0] / hh[1]));
}

// threshold: order and error set to 0 if error below threshold
template <class Field, class Range>
void PrintOrder(
    const Field& order, const Field& error, Range indices,
    Scal threshold = 1e-9) {
  Scal min_order = std::numeric_limits<Scal>::max();
  Scal max_order = -std::numeric_limits<Scal>::max();
  Scal min_error = 0; // error from min order location
  Scal max_error = 0; // error from max order location
  Scal mean_order = 0;
  Scal meanw = 0;
  Scal mean_error = 0;
  for (auto i : indices) {
    const Scal ord = Vect(order[i]).mean();
    if (ord < min_order) {
      min_order = ord;
      min_error = Vect(error[i]).mean();
    }
    if (ord > max_order) {
      max_order = ord;
      max_error = Vect(error[i]).mean();
    }
    mean_order += ord;
    mean_error += Vect(error[i]).mean();
    meanw += 1;
  }
  mean_order /= meanw;
  mean_error /= meanw;
  if (min_error < threshold) {
    min_order = 0;
    min_error = 0;
  }
  if (max_error < threshold) {
    max_order = 0;
    max_error = 0;
  }
  if (mean_error < threshold) {
    mean_order = 0;
    mean_error = 0;
  }
  printf(
      "order [error]: min=%6.3f [%10.3e]   max=%6.3f [%10.3e]   "
      "mean=%6.3f "
      "[%10.3e]\n",
      min_order, min_error, max_order, max_error, mean_order, mean_error);
}

// F derived from Func
template <class F, class Field>
void CalcAndPrintOrder(
    std::function<Field(const Func<typename F::Result>&, const M&)> estimator,
    std::function<Field(const Func<typename F::Result>&, const M&)> exact) {
  using Idx = typename Field::Idx;
  std::cout << "> func=" << F::GetName() << std::endl;
  const Scal h0 = 1e-3;
  const Field order = GetOrderField<F, Field, M>(estimator, exact);
  const Field error = GetErrorField<F, Field, M>(estimator, exact, h0);
  auto pm = CreateMesh(h0);
  auto& m = *pm;
  PrintOrder(order, error, m.GetRangeIn<Idx>());
}

// F derived from Func
template <class F, class T>
void CalcAndPrintOrder(
    std::function<FieldCell<T>(const Func<typename F::Result>&, const EB&)>
        estimator,
    std::function<FieldCell<T>(const Func<typename F::Result>&, const EB&)>
        exact) {
  using Field = FieldCell<T>;
  std::cout << "> func=" << F::GetName() << std::endl;
  const Scal h0 = 1e-3;
  const Field order = GetOrderField<F, Field, EB>(estimator, exact);
  const Field error = GetErrorField<F, Field, EB>(estimator, exact, h0);
  auto pm = CreateMesh(h0);
  auto peb = CreateEmbed(*pm);
  auto& m = *pm;
  auto& eb = *peb;
  printf(">> regular cells:\n   ");
  PrintOrder(order, error, MakeFilterIterator(m.Cells(), [&eb](IdxCell c) {
               return eb.GetType(c) == EB::Type::regular;
             }));
  printf(">> cut cells:\n   ");
  PrintOrder(order, error, MakeFilterIterator(m.Cells(), [&eb](IdxCell c) {
               return eb.GetType(c) == EB::Type::cut;
             }));
}

// F derived from Func
template <class F, class T>
void CalcAndPrintOrder(
    std::function<FieldEmbed<T>(const Func<typename F::Result>&, const EB&)>
        estimator,
    std::function<FieldEmbed<T>(const Func<typename F::Result>&, const EB&)>
        exact) {
  using Field = FieldEmbed<T>;
  std::cout << "> func=" << F::GetName() << std::endl;
  const Scal h0 = 1e-3;
  const Field order = GetOrderField<F, Field, EB>(estimator, exact);
  const Field error = GetErrorField<F, Field, EB>(estimator, exact, h0);
  auto pm = CreateMesh(h0);
  auto peb = CreateEmbed(*pm);
  auto& m = *pm;
  auto& eb = *peb;
  printf(">> regular faces between regular cells:\n   ");
  PrintOrder(order, error, MakeFilterIterator(m.Faces(), [&eb](IdxFace f) {
               return eb.GetType(eb.GetCell(f, 0)) == EB::Type::regular &&
                      eb.GetType(eb.GetCell(f, 1)) == EB::Type::regular;
             }));
  printf(">> regular faces near cut cell:\n   ");
  PrintOrder(order, error, MakeFilterIterator(m.Faces(), [&eb](IdxFace f) {
               return eb.GetType(f) == EB::Type::regular &&
                      (eb.GetType(eb.GetCell(f, 0)) == EB::Type::cut ||
                       eb.GetType(eb.GetCell(f, 1)) == EB::Type::cut);
             }));
  printf(">> cut faces:\n   ");
  PrintOrder(order, error, MakeFilterIterator(m.Faces(), [&eb](IdxFace f) {
               return eb.GetType(f) == EB::Type::cut;
             }));
  printf(">> embed faces:\n   ");
  PrintOrder(order, error, MakeFilterIterator(m.Cells(), [&eb](IdxCell c) {
               return eb.GetType(c) == EB::Type::cut;
             }));
}

void DumpField(const FieldCell<Scal>& fc, std::string filename, const M& m) {
  Dump(fc, m.GetIndexCells(), m.GetInBlockCells(), filename);
}

void DumpField(const FieldCell<Scal>& fc, std::string filename) {
  DumpField(fc, filename, *CreateMesh(1));
}

template <class Field, class MEB>
void VaryFunc(
    std::function<Field(const Func<Scal>&, const MEB&)> estimator,
    std::function<Field(const Func<Scal>&, const MEB&)> exact) {
  CalcAndPrintOrder<Linear>(estimator, exact);
  CalcAndPrintOrder<Trilinear>(estimator, exact);
  CalcAndPrintOrder<QuadraticBilinear>(estimator, exact);
  CalcAndPrintOrder<Sine>(estimator, exact);
}

template <class Field, class MEB>
void VaryFuncVect(
    std::function<Field(const Func<Vect>&, const MEB&)> estimator,
    std::function<Field(const Func<Vect>&, const MEB&)> exact) {
  CalcAndPrintOrder<FuncVect<QuadraticBilinear>>(estimator, exact);
  CalcAndPrintOrder<FuncVect<Sine>>(estimator, exact);
}

template <int dummy_ = 0>
void TestMesh() {
  std::cout << "\n" << __func__ << std::endl;
  std::cout << "\nExplViscous returning FieldCell<Vect> " << std::endl;
  auto estimator = [](const Func<Vect>& func, const M& m) {
    FieldCell<Vect> fcr(m, Vect(0));
    const FieldFace<Scal> ff_mu(m, 1);
    const auto fcu = Eval<FieldCell<Vect>>(func(), m);
    UFluid<M>::AppendExplViscous(fcr, fcu, {}, ff_mu, m);
    return fcr;
  };
  auto exact = [](const Func<Vect>& u, const M& m) {
    FieldCell<Vect> fcr(m, Vect(0));
    for (auto c : m.AllCells()) {
      const auto x = m.GetCenter(c);
      Vect r;
      r[0] = u.Dxx(0, 0)(x)[0] + u.Dxx(0, 1)(x)[1] + u.Dxx(0, 2)(x)[2];
      r[1] = u.Dxx(1, 0)(x)[0] + u.Dxx(1, 1)(x)[1] + u.Dxx(1, 2)(x)[2];
      r[2] = u.Dxx(2, 0)(x)[0] + u.Dxx(2, 1)(x)[1] + u.Dxx(2, 2)(x)[2];
      fcr[c] = r;
    }
    return fcr;
  };
  VaryFuncVect<FieldCell<Vect>, M>(estimator, exact);
}

template <int dummy_ = 0>
void DumpEmbedCsv() {
  auto pm = CreateMesh(1. / 16);
  auto peb = ConvertMesh<EB>(pm);
  auto& eb = *peb;
  std::ofstream out("eb_faces.csv");
  out << "x,y,z,nx,ny,nz,face,area,type\n";
  for (auto c : eb.CFaces()) {
    auto x = eb.GetFaceCenter(c);
    auto n = eb.GetNormal(c);
    out << x[0] << "," << x[1] << "," << x[2];
    out << "," << n[0] << "," << n[1] << "," << n[2];
    out << "," << 0;
    out << "," << eb.GetArea(c);
    out << "," << size_t(eb.GetType(c));
    out << "\n";
  }
  for (auto f : eb.Faces()) {
    auto x = eb.GetFaceCenter(f);
    auto n = eb.GetNormal(f);
    out << x[0] << "," << x[1] << "," << x[2];
    out << "," << n[0] << "," << n[1] << "," << n[2];
    out << "," << 1;
    out << "," << eb.GetArea(f);
    out << "," << size_t(eb.GetType(f));
    out << "\n";
  }
}

template <int dummy_ = 0>
void DumpEmbedPoly() {
  auto pm = CreateMesh(1. / 16);
  auto peb = ConvertMesh<EB>(pm);
  do {
    peb->DumpPoly();
  } while (pm->Pending());
}

template <int dummy_ = 0>
void DumpEmbedField() {
  auto pm = CreateMesh(1e-3);
  auto peb = ConvertMesh<EB>(pm);
  DumpField(Eval<FieldCell<Scal>>(Sine(0)(), *peb), "error_eb.dat", *pm);
}

template <int dummy_ = 0>
void TestEmbed() {
  std::cout << "\n" << __func__ << std::endl;

  if (1) {
    std::cout << "\n" << __func__ << std::endl;
    auto estimator = [](const Func<Scal>& func, const EB& eb) {
      auto fc = Eval<FieldCell<Scal>>(func(), eb);
      return UEB::Interpolate(fc, {}, eb);
    };
    auto exact = [](const Func<Scal>& func, const EB& eb) {
      return Eval<FieldEmbed<Scal>>(func(), eb);
    };
    std::cout << "\n"
              << "eb.Interpolate() returning FieldEmbed<Scal> " << std::endl;
    VaryFunc<FieldEmbed<Scal>, EB>(estimator, exact);
  }

  if (1) {
    auto estimator = [](const Func<Scal>& func, const EB& eb) {
      const auto fcu = Eval<FieldCell<Scal>>(func(), eb);
      return UEB::Gradient(fcu, {}, eb);
    };
    auto exact = [](const Func<Scal>& func, const EB& eb) {
      FieldEmbed<Scal> fe(eb, 0);
      for (auto c : eb.CFaces()) {
        const auto x = eb.GetFaceCenter(c);
        fe[c] = Gradient(func)(x).dot(eb.GetNormal(c));
      }
      for (auto f : eb.Faces()) {
        const auto x = eb.GetFaceCenter(f);
        fe[f] = Gradient(func)(x).dot(eb.GetNormal(f));
      }
      return fe;
    };
    std::cout << "\n"
              << "eb.GradientLimited() returning FieldEmbed<Scal>" << std::endl;
    VaryFunc<FieldEmbed<Scal>, EB>(estimator, exact);
  }

  if (1) {
    auto estimator = [](const Func<Scal>& func, const EB& eb) {
      auto fe = Eval<FieldEmbed<Scal>>(func(), eb);
      return UEB::GradientGauss(fe, eb);
    };
    auto estimator_lin = [](const Func<Scal>& func, const EB& eb) {
      auto fe = Eval<FieldEmbed<Scal>>(func(), eb);
      return UEB::GradientLinearFit(fe, eb);
    };
    auto exact = [](const Func<Scal>& func, const EB& eb) {
      return Eval<FieldCell<Vect>>(Gradient(func), eb);
    };
    std::cout << "\n"
              << "eb.GradientGauss() returning FieldCell<Vect>" << std::endl;
    VaryFunc<FieldCell<Vect>, EB>(estimator, exact);
    std::cout << "\n"
              << "eb.GradientLinearFit() returning FieldCell<Vect>"
              << std::endl;
    VaryFunc<FieldCell<Vect>, EB>(estimator_lin, exact);
  }

  if (1) {
    std::cout << "\n"
              << "eb.GradientGauss() returning FieldCell<Vect>" << std::endl;
    VaryFunc<FieldCell<Vect>, EB>(
        [](const Func<Scal>& func, const EB& eb) {
          auto fe = Eval<FieldEmbed<Scal>>(func(), eb);
          return UEB::GradientGauss(fe, eb);
        },
        [](const Func<Scal>& func, const EB& eb) {
          return Eval<FieldCell<Vect>>(Gradient(func), eb);
        });
  }

  using F = Quadratic;
  DumpField(
      GetOrderField<F, FieldCell<Scal>, EB>(
          [](const Func<Scal>& func, const EB& eb) {
            auto fe = Eval<FieldEmbed<Scal>>(func(), eb);
            return GetComponent(UEB::GradientGauss(fe, eb), 0);
          },
          [](const Func<Scal>& func, const EB& eb) {
            return Eval<FieldCell<Scal>>(func.Dx(0), eb);
          }),
      "order.dat");

  DumpField(
      Eval<F, FieldCell<Scal>, EB>(
          [](const Func<Scal>& func, const EB& eb) {
            auto fe = Eval<FieldEmbed<Scal>>(func(), eb);
            return GetComponent(UEB::GradientGauss(fe, eb), 0);
          },
          1e-3),
      "estimator.dat");

  DumpField(
      Eval<F, FieldCell<Scal>, EB>(
          [](const Func<Scal>& func, const EB& eb) {
            return Eval<FieldCell<Scal>>(func.Dx(0), eb);
          },
          1e-3),
      "exact.dat");

  DumpField(
      Eval<F, FieldCell<Scal>, EB>(
          [](const Func<Scal>& func, const EB& eb) {
            auto fe = Eval<FieldEmbed<Scal>>(func(), eb);
            return Sub(
                GetComponent(UEB::GradientGauss(fe, eb), 0),
                Eval<FieldCell<Scal>>(func.Dx(0), eb), eb);
          },
          1e-3),
      "error.dat");

  if (0) { // FIXME large error in cut cells
    std::cout << "\n"
              << "ExplViscous returning FieldCell<Vect> " << std::endl;
    auto estimator = [](const Func<Vect>& func, const EB& eb) {
      FieldCell<Vect> fcr(eb, Vect(0));
      const FieldFace<Scal> ff_mu(eb, 1);
      MapEmbed<BCond<Vect>> mebc;
      for (auto c : eb.CFaces()) {
        mebc[c] =
            BCond<Vect>(BCondType::dirichlet, 0, func()(eb.GetFaceCenter(c)));
      }
      const auto fcu = Eval<FieldCell<Vect>>(func(), eb.GetMesh());
      UFluid<M>::AppendExplViscous(fcr, fcu, {}, ff_mu, eb);
      return fcr;
    };
    auto exact = [](const Func<Vect>& u, const EB& eb) {
      FieldCell<Vect> fcr(eb, Vect(0));
      for (auto c : eb.AllCells()) {
        const auto x = eb.GetCellCenter(c);
        Vect r;
        r[0] = u.Dxx(0, 0)(x)[0] + u.Dxx(0, 1)(x)[1] + u.Dxx(0, 2)(x)[2];
        r[1] = u.Dxx(1, 0)(x)[0] + u.Dxx(1, 1)(x)[1] + u.Dxx(1, 2)(x)[2];
        r[2] = u.Dxx(2, 0)(x)[0] + u.Dxx(2, 1)(x)[1] + u.Dxx(2, 2)(x)[2];
        fcr[c] = r;
      }
      return fcr;
    };
    VaryFuncVect<FieldCell<Vect>, EB>(estimator, exact);
  }
}

template <int dummy_ = 0>
void TestEmbedUpwind() {
  if (1) {
    std::cout << "\n" << __func__ << std::endl;
    auto fluxdir = [](const EB& eb) {
      FieldEmbed<Scal> fev(eb); // flux direction, depends only on indices
      for (auto f : eb.Faces()) {
        fev[f] = std::sin(size_t(f));
      }
      for (auto c : eb.CFaces()) {
        fev[c] = std::sin(size_t(c));
      }
      return fev;
    };
    ConvSc sc = ConvSc::fou; // interpolation scheme
    auto estimator_sc = [&fluxdir, &sc](const Func<Scal>& func, const EB& eb) {
      MapCell<Scal> bc;
      for (auto c : eb.CFaces()) {
        bc[c] = func()(eb.GetFaceCenter(c));
      }
      const auto fcu = Eval<FieldCell<Scal>>(func(), eb.GetMesh());
      const auto feu = UEB::Interpolate(fcu, {}, eb);
      const auto fcg = UEB::GradientLinearFit(feu, eb);
      return UEB::InterpolateUpwind(fcu, {}, sc, fcg, fluxdir(eb), eb);
    };
    auto exact = [](const Func<Scal>& func, const EB& eb) {
      return Eval<FieldEmbed<Scal>>(func(), eb);
    };
    sc = ConvSc::sou;
    std::cout << "\n"
              << "eb.InterpolateUpwind() sc=sou returning FieldEmbed<Scal>"
              << std::endl;
    VaryFunc<FieldEmbed<Scal>, EB>(estimator_sc, exact);

    sc = ConvSc::quick;
    std::cout << "\n"
              << "eb.InterpolateUpwind() sc=quick returning FieldEmbed<Scal>"
              << std::endl;
    VaryFunc<FieldEmbed<Scal>, EB>(estimator_sc, exact);
  }
}

int main() {
  TestMesh();
  TestEmbed();
  TestEmbedUpwind();
  DumpEmbedPoly();
  DumpEmbedCsv();
}
