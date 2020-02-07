// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

template <class T>
T sqr(T a) {
  return a * a;
}

template <class Scal>
inline Scal GetNan() {
  return std::numeric_limits<Scal>::quiet_NaN();
}

template <>
inline double GetNan<double>() {
  return std::numeric_limits<double>::quiet_NaN();
}

template <class Scal>
bool IsFinite(Scal a) {
  return std::isfinite(a);
}

template <class Scal>
bool IsNan(Scal a) {
  return std::isnan(a);
}

namespace generic {

template <class Scal, size_t dim_>
class Vect {
 public:
  static constexpr size_t dim = dim_;
  using value_type = Scal;

 private:
  std::array<Scal, dim> comp_;

 public:
  friend void swap(Vect& first, Vect& second) {
    using std::swap;
    swap(first.comp_, second.comp_);
  }
  Vect() {}
  Vect(const Vect& vect) : comp_(vect.comp_) {}
  size_t size() const {
    return comp_.size();
  }
  explicit Vect(Scal value) {
    for (auto& a : comp_) {
      a = value;
    }
  }
  template <class... Args>
  explicit Vect(Scal first, Args... args) : comp_{{first, args...}} {
    constexpr size_t num_args = 1 + sizeof...(args);
    static_assert(
        num_args == dim,
        "Vect braced initializer must contain exactly 'dim' arguments");
  }
  template <class OtherScal>
  explicit Vect(const Vect<OtherScal, dim>& other) {
    for (size_t i = 0; i < dim; ++i) {
      comp_[i] = static_cast<Scal>(other[i]);
    }
  }
  template <class OtherScal>
  explicit Vect(const std::vector<OtherScal>& v) {
    // TODO dim instead of dim_ causes linker error
    for (size_t i = 0; i < dim; ++i) {
      comp_[i] = (i < v.size() ? static_cast<Scal>(v[i]) : 0);
    }
  }
  template <class OtherScal>
  explicit Vect(const std::array<OtherScal, dim>& v) {
    // TODO dim instead of dim_ causes linker error
    for (size_t i = 0; i < dim; ++i) {
      comp_[i] = static_cast<Scal>(v[i]);
    }
  }
  template <class OtherScal>
  explicit Vect(const OtherScal* v) {
    // TODO dim instead of dim_ causes linker error
    for (size_t i = 0; i < dim_; ++i) {
      comp_[i] = static_cast<Scal>(v[i]);
    }
  }
  Vect& operator=(Vect other) {
    comp_ = other.comp_;
    return *this;
  }
  Scal& operator[](size_t i) {
    return comp_[i];
  }
  const Scal& operator[](size_t i) const {
    return comp_[i];
  }
  Vect& operator+=(const Vect& vect) {
    for (size_t i = 0; i < dim; ++i) {
      comp_[i] += vect.comp_[i];
    }
    return *this;
  }
  Vect& operator-=(const Vect& other) {
    for (size_t i = 0; i < dim; ++i) {
      comp_[i] -= other.comp_[i];
    }
    return *this;
  }
  Vect& operator*=(Scal k) {
    for (size_t i = 0; i < dim; ++i) {
      comp_[i] *= k;
    }
    return *this;
  }
  Vect& operator/=(Scal k) {
    for (size_t i = 0; i < dim; ++i) {
      comp_[i] /= k;
    }
    return *this;
  }
  Vect operator+(Vect other) const {
    other += *this;
    return other;
  }
  Vect operator-(Vect other) const {
    Vect tmp(*this);
    tmp -= other;
    return tmp;
  }
  Vect operator-() const {
    Vect tmp(*this);
    for (size_t i = 0; i < dim; ++i) {
      tmp[i] = -tmp[i];
    }
    return tmp;
  }
  Vect operator*(Scal k) const {
    Vect tmp(*this);
    tmp *= k;
    return tmp;
  }
  Vect operator/(Scal k) const {
    Vect tmp(*this);
    tmp /= k;
    return tmp;
  }
  Vect& operator*=(const Vect& other) {
    for (size_t i = 0; i < dim; ++i) {
      comp_[i] *= other.comp_[i];
    }
    return *this;
  }
  Vect operator*(const Vect& other) const {
    Vect tmp(*this);
    tmp *= other;
    return tmp;
  }
  Vect& operator/=(const Vect& other) {
    for (size_t i = 0; i < dim; ++i) {
      comp_[i] /= other.comp_[i];
    }
    return *this;
  }
  Vect operator/(const Vect& other) const {
    Vect tmp(*this);
    tmp /= other;
    return tmp;
  }
  bool operator==(const Vect& other) const {
    for (size_t i = 0; i < dim; ++i) {
      if (!(comp_[i] == other.comp_[i])) {
        return false;
      }
    }
    return true;
  }
  bool operator!=(const Vect& other) const {
    return !(*this == other);
  }
  bool operator<(const Vect& other) const {
    for (size_t i = 0; i < dim; ++i) {
      if (!(comp_[i] < other.comp_[i])) {
        return false;
      }
    }
    return true;
  }
  bool operator<=(const Vect& other) const {
    for (size_t i = 0; i < dim; ++i) {
      if (!(comp_[i] <= other.comp_[i])) {
        return false;
      }
    }
    return true;
  }
  bool lexless(const Vect& o) const {
    return comp_ < o.comp_;
  }
  // TODO: remove, replace with Vect(0)
  static const Vect kZero;
  // TODO: remove, replace with Vect(1)
  static Vect GetUnit(size_t i) {
    Vect res = kZero;
    res[i] = 1;
    return res;
  }
  Scal sqrnorm() const {
    Scal res = 0;
    for (size_t i = 0; i < dim; ++i) {
      res += sqr(comp_[i]);
    }
    return res;
  }
  Scal norm() const {
    return std::sqrt(sqrnorm());
  }
  Scal dot(const Vect& other) const {
    Scal sum = 0;
    for (size_t i = 0; i < dim; ++i) {
      sum += comp_[i] * other.comp_[i];
    }
    return sum;
  }
  Scal cross_third(const Vect& other) const {
    return comp_[0] * other.comp_[1] - comp_[1] * other.comp_[0];
  }
  Vect cross(const Vect& other) const {
    const Vect& a = *this;
    const Vect& b = other;
    return Vect(
        a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]);
  }
  Scal dist(Vect other) const {
    other -= *this;
    return other.norm();
  }
  Scal sqrdist(Vect o) const {
    o -= *this;
    return o.sqrnorm();
  }
  Scal sum() const {
    Scal r = 0;
    for (size_t i = 0; i < dim; ++i) {
      r += comp_[i];
    }
    return r;
  }
  Scal mean() const {
    return sum() / dim;
  }
  Scal prod() const {
    Scal r = comp_[0];
    for (size_t i = 1; i < dim; ++i) {
      r *= comp_[i];
    }
    return r;
  }
  Scal norm1() const {
    Scal r = 0;
    for (size_t i = 0; i < dim; ++i) {
      r += std::abs(comp_[i]);
    }
    return r;
  }
  Scal norminf() const {
    Scal r = 0;
    for (size_t i = 0; i < dim; ++i) {
      r = std::max(r, std::abs(comp_[i]));
    }
    return r;
  }
  Scal max() const {
    Scal r = comp_[0];
    for (size_t i = 1; i < dim; ++i) {
      r = std::max(r, comp_[i]);
    }
    return r;
  }
  Scal min() const {
    Scal r = comp_[0];
    for (size_t i = 1; i < dim; ++i) {
      r = std::min(r, comp_[i]);
    }
    return r;
  }
  size_t argmax() const {
    size_t r = 0;
    for (size_t i = 1; i < dim; ++i) {
      if (comp_[i] > comp_[r]) {
        r = i;
      }
    }
    return r;
  }
  size_t argmin() const {
    size_t r = 0;
    for (size_t i = 1; i < dim; ++i) {
      if (comp_[i] < comp_[r]) {
        r = i;
      }
    }
    return r;
  }
  Vect abs() const {
    Vect r = *this;
    for (size_t i = 0; i < dim; ++i) {
      r[i] = std::abs(r[i]);
    }
    return r;
  }
  Vect max(Vect o) const {
    for (size_t i = 0; i < dim; ++i) {
      o[i] = std::max(comp_[i], o[i]);
    }
    return o;
  }
  Vect min(Vect o) const {
    for (size_t i = 0; i < dim; ++i) {
      o[i] = std::min(comp_[i], o[i]);
    }
    return o;
  }
  Vect clip(const Vect& v0, const Vect& v1) const {
    return (*this).max(v0).min(v1);
  }
  // TODO: revise, may lead to undesired conversion
  template <class T = Scal>
  operator std::vector<T>() const {
    return std::vector<T>(comp_.begin(), comp_.end());
  }
  template <class T = Scal>
  operator std::array<T, dim>() const {
    std::array<T, dim> r;
    for (size_t i = 0; i < dim; ++i) {
      r[i] = static_cast<T>(comp_[i]);
    }
    return r;
  }
  class LexLess {
   public:
    bool operator()(Vect a, Vect b) const {
      return a.lexless(b);
    }
  };
  friend std::ostream& operator<<(std::ostream& out, const Vect& v) {
    out << "(";
    for (size_t i = 0; i < dim; ++i) {
      if (i != 0) {
        out << ",";
      }
      out << v[i];
    }
    out << ")";
    return out;
  }

  friend std::istream& operator>>(std::istream& in, Vect& v) {
    for (size_t i = 0; i < dim; ++i) {
      in >> v[i];
    }
    return in;
  }
};

template <class Scal, size_t dim>
const Vect<Scal, dim> Vect<Scal, dim>::kZero =
    Vect<Scal, dim>(static_cast<Scal>(0.));

} // namespace generic

template <class Vect>
struct Rect {
  static constexpr size_t dim = Vect::dim;

  Vect lb, rt;
  Rect() {}
  Rect(const Vect& lb, const Vect& rt) : lb(lb), rt(rt) {}
  bool IsInside(Vect x) const {
    for (size_t i = 0; i < dim; ++i) {
      if (x[i] < lb[i] || rt[i] < x[i]) {
        return false;
      }
    }
    return true;
  }
  Vect GetDimensions() const {
    return rt - lb;
  }
};

template <>
inline generic::Vect<double, 3> GetNan<generic::Vect<double, 3>>() {
  return generic::Vect<double, 3>(GetNan<double>());
}
