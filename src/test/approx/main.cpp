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
#include "solver/cond.h"
#include "solver/embed.h"
#include "solver/solver.h"

using M = MeshStructured<double, 3>;
using Scal = typename M::Scal;
using Vect = typename M::Vect;
using MIdx = typename M::MIdx;
auto Dirs = GRange<size_t>(3);
using EB = Embed<M>;

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
  std::default_random_engine gen_;
};

class Quadratic : public Func<Scal> {
 public:
  explicit Quadratic(int seed)
      : Func<Scal>(seed), origin_(Vect(1) + Vect(G(), G(), G()) * 0.1) {}
  virtual ~Quadratic() = default;
  std::function<Scal(Vect)> operator()() const override {
    return [this](Vect v) {
      v += origin_;
      const Scal x(v[0]), y(v[1]), z(v[2]);
      return x * x + y * y + z * z;
    };
  }
  std::function<Scal(Vect)> Dx(size_t d) const override {
    return [d, this](Vect v) {
      v += origin_;
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

class Sine : public Func<Scal> {
 public:
  explicit Sine(int seed)
      : Func<Scal>(seed)
      , freq_(Vect(G(), G(), G()) + Vect(10))
      , phase_(Vect(0.5) + Vect(G(), G(), G()) * 10) {}
  virtual ~Sine() = default;
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
class FuncVect : public Func<Vect> {
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
    Field r(m);
    for (auto i : m.template GetAll<Idx>()) {
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
    for (auto i : m.template GetAll<Idx>()) {
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
    for (auto i : m.template GetAll<Idx>()) {
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
  return std::make_unique<M>(
      InitUniformMesh<M>(dom, MIdx(0), size, 2, true, true, size, 0));
}

// TODO: rename GetInBlockCells to GetBlockInCells()
// TODO: rename SuCells to Cells(1)
// TODO: rename AllCells to Cells(2)

std::unique_ptr<EB> CreateEmbed(M& m) {
  auto peb = std::make_unique<EB>(m, 0);
  FieldNode<Scal> fnl(m);
  auto block = m.GetInBlockCells().GetSize();
  auto h = m.GetCellSize();
  for (auto n : m.AllNodes()) {
    const auto x = m.GetNode(n) / h / Vect(block);
    fnl[n] = 0.4 - Vect(0.5).dist(x);
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

template <class T, class Idx, class MEB>
Scal Norm1(const GField<T, Idx>& u, const MEB& meb) {
  Scal sum = 0;
  Scal sumv = 0;
  for (auto i : GRange<Idx>(meb)) {
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
        return std::log(e0 / e1) / std::log(h0 / h1);
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

// F derived from Func
template <class F, class Field, class MEB>
void PrintOrder(
    std::function<Field(const Func<typename F::Result>&, const MEB&)> estimator,
    std::function<Field(const Func<typename F::Result>&, const MEB&)> exact) {
  std::cout << "func=" << F::GetName() << std::endl;
  auto t = CalcOrder<F, Field>(estimator, exact);
  printf(
      "order=%5.3f   coarse=%.5e   fine=%.5e\n", std::get<2>(t), std::get<0>(t),
      std::get<1>(t));
}

void DumpField(const FieldCell<Scal>& fc, std::string filename, const M& m) {
  Dump(fc, m.GetIndexCells(), m.GetInBlockCells(), filename);
}

template <class Field, class MEB>
void VaryFunc(
    std::function<Field(const Func<Scal>&, const MEB&)> estimator,
    std::function<Field(const Func<Scal>&, const MEB&)> exact) {
  PrintOrder<Quadratic, Field>(estimator, exact);
  PrintOrder<Sine, Field>(estimator, exact);
}

template <class Field, class MEB>
void VaryFuncVect(
    std::function<Field(const Func<Vect>&, const MEB&)> estimator,
    std::function<Field(const Func<Vect>&, const MEB&)> exact) {
  PrintOrder<FuncVect<Quadratic>, Field>(estimator, exact);
  PrintOrder<FuncVect<Sine>, Field>(estimator, exact);
}

void TestMesh() {
  std::cout << "\n"
            << "GradientI() returning FieldFace<Scal>" << std::endl;
  VaryFunc<FieldFace<Scal>, M>(
      [](const Func<Scal>& func, const M& m) {
        const auto fcu = Eval<FieldCell<Scal>>(func(), m);
        FieldFace<Scal> ffg(m);
        GradientI(fcu, m, ffg);
        return ffg;
      },
      [](const Func<Scal>& func, const M& m) {
        FieldFace<Scal> ffg(m, 0);
        for (auto f : m.AllFaces()) {
          const auto x = m.GetCenter(f);
          ffg[f] = Gradient(func)(x).dot(m.GetNormal(f));
        }
        return ffg;
      });

  std::cout << "\n"
            << "Average() returning FieldCell<Scal> " << std::endl;
  VaryFunc<FieldCell<Scal>, M>(
      [](const Func<Scal>& func, const M& m) {
        const auto ffu = Eval<FieldFace<Scal>>(func(), m);
        const auto fcu = Average(ffu, m);
        return fcu;
      },
      [](const Func<Scal>& func, const M& m) {
        return Eval<FieldCell<Scal>>(func(), m);
      });

  std::cout << "\n"
            << "Laplace returning FieldCell<Scal> " << std::endl;
  VaryFunc<FieldCell<Scal>, M>(
      [](const Func<Scal>& func, const M& m) {
        const auto fcu = Eval<FieldCell<Scal>>(func(), m);
        FieldFace<Scal> ffg(m);
        GradientI(fcu, m, ffg);
        FieldCell<Scal> fcl(m);
        for (auto c : m.Cells()) {
          Scal sum = 0;
          for (auto nci : m.Nci(c)) {
            const IdxFace f = m.GetFace(c, nci);
            sum += ffg[f] * m.GetArea(f) * m.GetOutwardFactor(c, nci);
          }
          fcl[c] = sum / m.GetVolume(c);
        }
        return fcl;
      },
      [](const Func<Scal>& func, const M& m) {
        FieldCell<Scal> r(m, 0);
        for (auto c : m.AllCells()) {
          const auto x = m.GetCenter(c);
          const Scal uxx = func.Dxx(0, 0)(x);
          const Scal uyy = func.Dxx(1, 1)(x);
          const Scal uzz = func.Dxx(2, 2)(x);
          r[c] = uxx + uyy + uzz;
        }
        return r;
      });

  std::cout << "\n"
            << "Average() returning FieldCell<Vect> " << std::endl;
  VaryFuncVect<FieldCell<Vect>, M>(
      [](const Func<Vect>& func, const M& m) {
        const auto ffu = Eval<FieldFace<Vect>>(func(), m);
        const auto fcu = Average(ffu, m);
        return fcu;
      },
      [](const Func<Vect>& func, const M& m) {
        return Eval<FieldCell<Vect>>(func(), m);
      });

  std::cout << "\n"
            << "Laplace returning FieldCell<Vect> " << std::endl;
  VaryFuncVect<FieldCell<Vect>, M>(
      [](const Func<Vect>& func, const M& m) {
        const auto fcu = Eval<FieldCell<Vect>>(func(), m);
        FieldFace<Vect> ffg(m);
        GradientI(fcu, m, ffg);
        FieldCell<Vect> fcl(m);
        for (auto c : m.Cells()) {
          Vect sum(0);
          for (auto nci : m.Nci(c)) {
            const IdxFace f = m.GetFace(c, nci);
            sum += ffg[f] * m.GetArea(f) * m.GetOutwardFactor(c, nci);
          }
          fcl[c] = sum / m.GetVolume(c);
        }
        return fcl;
      },
      [](const Func<Vect>& func, const M& m) {
        FieldCell<Vect> r(m, Vect(0));
        for (auto c : m.AllCells()) {
          const auto x = m.GetCenter(c);
          const Vect uxx = func.Dxx(0, 0)(x);
          const Vect uyy = func.Dxx(1, 1)(x);
          const Vect uzz = func.Dxx(2, 2)(x);
          r[c] = uxx + uyy + uzz;
        }
        return r;
      });

  std::cout << "\n"
            << "ExplViscous returning FieldCell<Vect> " << std::endl;
  VaryFuncVect<FieldCell<Vect>, M>(
      [](const Func<Vect>& func, const M& m) {
        FieldCell<Vect> fcr(m, Vect(0));
        const auto fcu = Eval<FieldCell<Vect>>(func(), m);
        auto ffu = Interpolate(fcu, MapCondFace(), m);
        for (auto d : GRange<size_t>(Vect::dim)) {
          auto fcg = Gradient(GetComponent(ffu, d), m);
          auto ffg = Interpolate(fcg, MapCondFace(), m);
          for (auto c : m.Cells()) {
            Vect s(0);
            for (auto q : m.Nci(c)) {
              const IdxFace f = m.GetFace(c, q);
              s += ffg[f] * m.GetOutwardSurface(c, q)[d];
            }
            fcr[c] += s / m.GetVolume(c);
          }
        }

        return fcr;
      },
      [](const Func<Vect>& u, const M& m) {
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
      });

  std::cout << "\n"
            << "Laplace dump error field FieldScal<Scal>" << std::endl;
  DumpField(
      GetErrorField<Sine, FieldCell<Scal>, M>(
          [](const Func<Scal>& func, const M& m) {
            const auto fcu = Eval<FieldCell<Scal>>(func(), m);
            FieldFace<Scal> ffg(m);
            GradientI(fcu, m, ffg);
            FieldCell<Scal> fcl(m);
            for (auto c : m.Cells()) {
              Scal sum = 0;
              for (auto nci : m.Nci(c)) {
                const IdxFace f = m.GetFace(c, nci);
                sum += ffg[f] * m.GetArea(f) * m.GetOutwardFactor(c, nci);
              }
              fcl[c] = sum / m.GetVolume(c);
            }
            return fcl;
          },
          [](const Func<Scal>& func, const M& m) {
            FieldCell<Scal> r(m, 0);
            for (auto c : m.AllCells()) {
              const auto x = m.GetCenter(c);
              const Scal uxx = func.Dxx(0, 0)(x);
              const Scal uyy = func.Dxx(1, 1)(x);
              const Scal uzz = func.Dxx(2, 2)(x);
              r[c] = uxx + uyy + uzz;
            }
            return r;
          },
          1e-3),
      "error.dat", *CreateMesh(1));
}

void DumpEmbedPoly() {
  auto pm = CreateMesh(1. / 16);
  auto peb = ConvertMesh<EB>(pm);
  do {
    peb->DumpPoly();
  } while (pm->Pending());
}

void DumpEmbedField() {
  auto pm = CreateMesh(1e-3);
  auto peb = ConvertMesh<EB>(pm);
  DumpField(
      Eval<FieldCell<Scal>>(
          //[](Vect x) { return x[0]; },
          Sine(0)(), *peb),
      "error_eb.dat", *pm);
}

void TestEmbed() {
  DumpEmbedPoly();
  std::cout << "\n"
            << "eb.Interpolate() returning FieldEmbed<Scal> " << std::endl;
  VaryFunc<FieldEmbed<Scal>, EB>(
      [](const Func<Scal>& func, const EB& eb) {
        auto fc = Eval<FieldCell<Scal>>(func(), eb);
        return eb.Interpolate(fc, MapCondFace(), 1, 0.);
      },
      [](const Func<Scal>& func, const EB& eb) {
        return Eval<FieldEmbed<Scal>>(func(), eb);
      });

  std::cout << "\n"
            << "eb.Gradient() returning FieldEmbed<Scal>" << std::endl;
  VaryFunc<FieldEmbed<Scal>, EB>(
      [](const Func<Scal>& func, const EB& eb) {
        const auto fcu = Eval<FieldCell<Scal>>(func(), eb);
        return eb.Gradient(fcu, MapCondFace(), 1, 0.);
      },
      [](const Func<Scal>& func, const EB& eb) {
        FieldEmbed<Scal> fe(eb, 0);
        for (auto c : eb.CFaces()) {
          fe[c] = 0;
        }
        for (auto f : eb.Faces()) {
          const auto x = eb.GetFaceCenter(f);
          fe[f] = Gradient(func)(x).dot(eb.GetNormal(f));
        }
        return fe;
      });

  std::cout << "\n"
            << "eb.Gradient() returning FieldCell<Vect>" << std::endl;
  VaryFunc<FieldCell<Vect>, EB>(
      [](const Func<Scal>& func, const EB& eb) {
        auto fe = Eval<FieldEmbed<Scal>>(func(), eb);
        return eb.Gradient(fe);
      },
      [](const Func<Scal>& func, const EB& eb) {
        return Eval<FieldCell<Vect>>(Gradient(func), eb);
      });

  DumpField(
      GetOrderField<Sine, FieldCell<Scal>, EB>(
          [](const Func<Scal>& func, const EB& eb) {
            auto fe = Eval<FieldEmbed<Scal>>(func(), eb);
            return GetComponent(eb.Gradient(fe), 0);
          },
          [](const Func<Scal>& func, const EB& eb) {
            return Eval<FieldCell<Scal>>(func.Dx(0), eb);
          }),
      "error_eb.dat", *CreateMesh(1));
}

int main() {
  // TestMesh();
  TestEmbed();
}
