/*
 * vect.hpp
 *
 *  Created on: Jan 25, 2016
 *      Author: Petr Karnakov (petr.karnakov@gmail.com)
 */

#pragma once

#include <array>
#include <initializer_list>
#include <iostream>
#include <cmath>
#include <cassert>

namespace geom {

template <class T>
T sqr(T a) {
  return a * a;
}

template <class Scal, size_t dimarg>
class Vect {
 public:
  static constexpr size_t dim = dimarg;
  using value_type = Scal;

 private:
  std::array<Scal, dim> comp_;

 public:
  friend void swap(Vect& first, Vect& second) {
    using std::swap;
    swap(first.comp_, second.comp_);
  }
  Vect() {}
  Vect(const Vect& vect)
      : comp_(vect.comp_)
  {}
  size_t size() const {
    return comp_.size();
  }
  explicit Vect(Scal value) {
    std::fill(comp_.begin(), comp_.end(), value);
  }
  template <class... Args>
  explicit Vect(Scal first, Args... args)
      : comp_{{first, args...}}
  {
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
  Vect& operator=(Vect other) {
    comp_ = other.comp_;
    return *this;
  }
  Scal& operator[](size_t i) {
#ifdef __RANGE_CHECK
    assert(i >=0 && i < comp_.size());
#endif
    return comp_[i];
  }
  const Scal& operator[](size_t i) const {
#ifdef __RANGE_CHECK
    assert(i >=0 && i < comp_.size());
#endif
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
  static const Vect kZero;
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
    return Vect(a[1]*b[2]-a[2]*b[1],
                a[2]*b[0]-a[0]*b[2],
                a[0]*b[1]-a[1]*b[0]);
  }
  Scal dist(Vect other) const {
    other -= *this;
    return other.norm();
  }
};

template <class Scal, size_t dim>
const Vect<Scal, dim> Vect<Scal, dim>::kZero =
    Vect<Scal, dim>(static_cast<Scal>(0.));

template <class Scal, size_t dim>
std::ostream& operator<<(std::ostream& out, const Vect<Scal, dim>& vect) {
  out << "(";
  for (size_t i = 0; i < dim; ++i) {
    if (i != 0) {
      out << ",";
    }
    out << vect[i];
  }
  out << ")";
  return out;
}

template <class Scal, size_t dim>
std::istream& operator>>(std::istream& in, Vect<Scal, dim>& vect) {
  for (size_t i = 0; i < dim; ++i) {
    in >> vect[i];
  }
  return in;
}

template <class Vect>
struct Rect {
  size_t dim = Vect::dim;
  Vect lb, rt;
  Rect() {}
  Rect(const Vect& lb, const Vect& rt)
      : lb(lb), rt(rt)
  {}
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

} // namespace geom
