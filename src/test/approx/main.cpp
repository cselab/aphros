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

template <class T>
class Func {
 public:
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

template <class T, class Idx>
GField<T, Idx> Add(
    const GField<T, Idx>& u0, const GField<T, Idx>& u1, const M& m) {
  GField<T, Idx> d(m);
  for (auto i : m.template GetAll<Idx>()) {
    d[i] = u0[i] + u1[i];
  }
  return d;
}

template <class T, class Idx>
GField<T, Idx> Sub(
    const GField<T, Idx>& u0, const GField<T, Idx>& u1, const M& m) {
  GField<T, Idx> r(m);
  for (auto i : m.template GetAll<Idx>()) {
    r[i] = u0[i] - u1[i];
  }
  return r;
}

template <class T, class Idx>
GField<T, Idx> Mul(const GField<T, Idx>& u, Scal k, const M& m) {
  GField<T, Idx> r(m);
  for (auto i : m.template GetAll<Idx>()) {
    r[i] = u[i] * k;
  }
  return r;
}

Embed<M> GetEmbed(M& m) {
  Embed<M> eb(m);
  FieldNode<Scal> fnl(m);
  for (auto n : m.AllNodes()) {
    const auto x = m.GetNode(n);
    fnl[n] = (0.4 - (Vect(0.5, 0.5, 0.) - Vect(x[0], x[1], 0.)).norm());
  }
  eb.Init(fnl);
  return eb;
}

template <class T>
FieldEmbed<T> GetFieldEmbed(
    const std::function<T(Vect)>& u, const Embed<M>& eb) {
  FieldEmbed<T> fe(eb.GetMesh());
  for (auto f : eb.Faces()) {
    fe[f] = u(eb.GetFaceCenter(f));
  }
  for (auto c : eb.CFaces()) {
    fe[c] = u(eb.GetFaceCenter(c));
  }
  return fe;
}

template <class T>
FieldCell<T> GetFieldCell(const std::function<T(Vect)>& u, const Embed<M>& eb) {
  FieldCell<T> fc(eb.GetMesh());
  for (auto c : eb.Cells()) {
    fc[c] = u(eb.GetCellCenter(c));
  }
  return fc;
}

template <class T>
FieldEmbed<T> Add(
    const FieldEmbed<T>& u0, const FieldEmbed<T>& u1, const Embed<M>& eb) {
  FieldEmbed<T> r(eb.GetMesh());
  for (auto c : eb.CFaces()) {
    r[c] = u0[c] + u1[c];
  }
  for (auto f : eb.Faces()) {
    r[f] = u0[f] + u1[f];
  }
  return r;
}

template <class T>
FieldEmbed<T> Sub(
    const FieldEmbed<T>& u0, const FieldEmbed<T>& u1, const Embed<M>& eb) {
  FieldEmbed<T> r(eb.GetMesh());
  for (auto c : eb.CFaces()) {
    r[c] = u0[c] - u1[c];
  }
  for (auto f : eb.Faces()) {
    r[f] = u0[f] - u1[f];
  }
  return r;
}

template <class T>
FieldEmbed<T> Mul(const FieldEmbed<T>& u, Scal k, const Embed<M>& eb) {
  FieldEmbed<T> r(eb.GetMesh());
  for (auto c : eb.CFaces()) {
    r[c] = u[c] * k;
  }
  for (auto f : eb.Faces()) {
    r[f] = u[f] * k;
  }
  return r;
}

Scal Abs(Scal a) {
  return std::abs(a);
}

Vect Abs(Vect a) {
  return a.abs();
}

Scal Norm1(Scal a) {
  return std::abs(a);
}

Scal Norm1(Vect a) {
  return a.norm1();
}

template <class T, class Idx>
GField<T, Idx> Abs(const GField<T, Idx>& u, const M& m) {
  GField<T, Idx> r(m);
  for (auto i : m.template GetAll<Idx>()) {
    r[i] = Abs(u[i]);
  }
  return r;
}

template <class T, class Idx, class M>
Scal Norm1(const GField<T, Idx>& u, const M& m) {
  Scal sum = 0;
  Scal sumv = 0;
  for (auto i : m.template GetIn<Idx>()) {
    sum += Norm1(u[i]);
    sumv += 1;
  }
  return sum / sumv;
}

template <class T>
FieldEmbed<T> Abs(const FieldEmbed<T>& u, const Embed<M>& eb) {
  FieldEmbed<T> r(eb.GetMesh());
  for (auto c : eb.CFaces()) {
    r[c] = Abs(u[c]);
  }
  for (auto f : eb.Faces()) {
    r[f] = Abs(u[f]);
  }
  return r;
}

template <class T>
Scal Norm1(const FieldEmbed<T>& u, const Embed<M>& eb) {
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
// F: function to evaluate, derived from Func<FT>, constructable from int
// func: function to evaluate
// estimator: estimator of field
// exact: exact field
// h: mesh step
// Returns:
// tuple(error0, error1, order)
template <class F, class FT, class Field>
Field GetErrorField(
    std::function<Field(const Func<FT>&, const M&)> estimator,
    std::function<Field(const Func<FT>&, const M&)> exact, Scal h) {
  const size_t nsamp = 10;
  const auto m = GetMesh(h);
  using Value = typename Field::Value;
  Field avg(m, Value(0));
  for (size_t seed = 0; seed < nsamp; ++seed) {
    const F func(seed);
    auto f0 = estimator(func, m);
    auto f1 = exact(func, m);
    avg = Add(avg, Abs(Sub(f0, f1, m), m), m);
  }
  avg = Mul(avg, 1. / nsamp, m);
  return avg;
}

template <class F, class FT, class Field>
Field GetErrorField(
    std::function<Field(const Func<FT>&, const Embed<M>&)> estimator,
    std::function<Field(const Func<FT>&, const Embed<M>&)> exact, Scal h) {
  const size_t nsamp = 10;
  auto m = GetMesh(h);
  const auto eb = GetEmbed(m);
  using Value = typename Field::Value;
  Field avg(m, Value(0));
  for (size_t seed = 0; seed < nsamp; ++seed) {
    const F func(seed);
    auto f0 = estimator(func, eb);
    auto f1 = exact(func, eb);
    avg = Add(avg, Abs(Sub(f0, f1, eb), eb), eb);
  }
  avg = Mul(avg, 1. / nsamp, eb);
  return avg;
}

/// Computes mean error field of estimator compared to exact values.
/// @param F function to evaluate, derived from Func<FT>, constructable from int
/// @param func function to evaluate
/// @param estimator estimator of field
/// @param exact exact field
/// @param h mesh step
/// @return tuple(error0, error1, order)

// Evaluates error on two meshes and estimates the order of accuracy.
// F: function to evaluate, derived from Func<FT>, constructable from int
// func: function to evaluate
// estimator: estimator of field
// exact: exact field
// Returns:
// tuple(error0, error1, order)
template <class F, class FT, class Field>
std::tuple<Scal, Scal, Scal> CalcOrder(
    std::function<Field(const Func<FT>&, const M&)> estimator,
    std::function<Field(const Func<FT>&, const M&)> exact) {
  const Scal h0 = 1e-3;
  const std::vector<Scal> hh = {h0, h0 * 0.5};
  std::vector<Scal> ee(hh.size());
  for (size_t i = 0; i < hh.size(); ++i) {
    ee[i] = Norm1(
        GetErrorField<F, FT, Field>(estimator, exact, hh[i]), GetMesh(hh[i]));
  }
  return std::make_tuple(
      ee[0], ee[1], std::log(ee[0] / ee[1]) / std::log(hh[0] / hh[1]));
}

// F derived from Func<FT>
template <class F, class FT, class Field>
void PrintOrder(
    std::function<Field(const Func<FT>&, const M&)> estimator,
    std::function<Field(const Func<FT>&, const M&)> exact) {
  std::cout << "func=" << F::GetName() << std::endl;
  auto t = CalcOrder<F, FT, Field>(estimator, exact);
  printf(
      "order=%5.3f   coarse=%.5e   fine=%.5e\n", std::get<2>(t), std::get<0>(t),
      std::get<1>(t));
}

void DumpField(const FieldCell<Scal>& fc, std::string filename, const M& m) {
  Dump(fc, m.GetIndexCells(), m.GetInBlockCells(), filename);
}

template <class Field>
void VaryFunc(
    std::function<Field(const Func<Scal>&, const M&)> estimator,
    std::function<Field(const Func<Scal>&, const M&)> exact) {
  PrintOrder<Quadratic, Scal, Field>(estimator, exact);
  PrintOrder<Sine, Scal, Field>(estimator, exact);
}

template <class Field>
void VaryFuncVect(
    std::function<Field(const Func<Vect>&, const M&)> estimator,
    std::function<Field(const Func<Vect>&, const M&)> exact) {
  PrintOrder<FuncVect<Quadratic>, Vect, Field>(estimator, exact);
  PrintOrder<FuncVect<Sine>, Vect, Field>(estimator, exact);
}

void TestMesh() {
  std::cout << "\n"
            << "GradientI() returning FieldFace<Scal>" << std::endl;
  VaryFunc<FieldFace<Scal>>(
      [](const Func<Scal>& func, const M& m) {
        const auto fcu = GetField<Scal, IdxCell>(func(), m);
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
  VaryFunc<FieldCell<Scal>>(
      [](const Func<Scal>& func, const M& m) {
        const auto ffu = GetField<Scal, IdxFace>(func(), m);
        const auto fcu = Average(ffu, m);
        return fcu;
      },
      [](const Func<Scal>& func, const M& m) {
        return GetField<Scal, IdxCell>(func(), m);
      });

  std::cout << "\n"
            << "Laplace returning FieldCell<Scal> " << std::endl;
  VaryFunc<FieldCell<Scal>>(
      [](const Func<Scal>& func, const M& m) {
        const auto fcu = GetField<Scal, IdxCell>(func(), m);
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
  VaryFuncVect<FieldCell<Vect>>(
      [](const Func<Vect>& func, const M& m) {
        const auto ffu = GetField<Vect, IdxFace>(func(), m);
        const auto fcu = Average(ffu, m);
        return fcu;
      },
      [](const Func<Vect>& func, const M& m) {
        return GetField<Vect, IdxCell>(func(), m);
      });

  std::cout << "\n"
            << "Laplace returning FieldCell<Vect> " << std::endl;
  VaryFuncVect<FieldCell<Vect>>(
      [](const Func<Vect>& func, const M& m) {
        const auto fcu = GetField<Vect, IdxCell>(func(), m);
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
  VaryFuncVect<FieldCell<Vect>>(
      [](const Func<Vect>& func, const M& m) {
        FieldCell<Vect> fcr(m, Vect(0));
        const auto fcu = GetField<Vect, IdxCell>(func(), m);
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
      GetErrorField<Sine, Scal, FieldCell<Scal>>(
          [](const Func<Scal>& func, const M& m) {
            const auto fcu = GetField<Scal, IdxCell>(func(), m);
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
      "error.dat", GetMesh(1));
}

void TestEmbed() {
  GetErrorField<Sine, Scal, FieldEmbed<Scal>>(
      [](const Func<Scal>& func, const Embed<M>& eb) {
        FieldEmbed<Scal> feu(eb.GetMesh(), 0);
        return feu;
      },
      [](const Func<Scal>& func, const Embed<M>& eb) {
        FieldEmbed<Scal> feu(eb.GetMesh(), 1);
        return feu;
      },
      1e-3);
  /*
  std::cout << "\n"
            << "eb.Interpolate() returning FieldCell<Scal>" << std::endl;
  VaryFunc<Scal, IdxCell>(
      [](const Func<Scal>& func, const M& m) {
        const auto fcu = GetField<Scal, IdxCell>(func(), m);
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
      */
}

int main() {
  //TestMesh();
  TestEmbed();
}
