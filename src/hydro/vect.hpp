#pragma once

#include <array>
#include <initializer_list>
#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>



template <class T>
T sqr(T a) {
  return a * a;
}

// TODO: revise __RANGE_CHECK

template <class Scal, size_t dimarg>
class GVect {
 public:
  static constexpr size_t dim = dimarg;
  using value_type = Scal;

 private:
  std::array<Scal, dim> comp_;

 public:
  friend void swap(GVect& first, GVect& second) {
    using std::swap;
    swap(first.comp_, second.comp_);
  }
  GVect() {}
  GVect(const GVect& vect)
      : comp_(vect.comp_)
  {}
  size_t size() const {
    return comp_.size();
  }
  explicit GVect(Scal value) {
    std::fill(comp_.begin(), comp_.end(), value);
  }
  template <class... Args>
  explicit GVect(Scal first, Args... args)
      : comp_{{first, args...}}
  {
    constexpr size_t num_args = 1 + sizeof...(args);
    static_assert(
        num_args == dim,
        "Vect braced initializer must contain exactly 'dim' arguments");
  }
  template <class OtherScal>
  explicit GVect(const GVect<OtherScal, dim>& other) {
    for (size_t i = 0; i < dim; ++i) {
      comp_[i] = static_cast<Scal>(other[i]);
    }
  }
  template <class OtherScal>
  explicit GVect(const std::vector<OtherScal>& v) {
    // TODO dim instead of dimarg causes linker error
    for (size_t i = 0; i < std::min<size_t>(dimarg, v.size()); ++i) {
      comp_[i] = static_cast<Scal>(v[i]);
    }
  }
  template <class OtherScal>
  explicit GVect(const std::array<OtherScal, dim>& v) {
    // TODO dim instead of dimarg causes linker error
    for (size_t i = 0; i < std::min<size_t>(dimarg, v.size()); ++i) {
      comp_[i] = static_cast<Scal>(v[i]);
    }
  }
  template <class OtherScal>
  explicit GVect(const OtherScal* v) {
    // TODO dim instead of dimarg causes linker error
    for (size_t i = 0; i < dimarg; ++i) {
      comp_[i] = static_cast<Scal>(v[i]);
    }
  }
  GVect& operator=(GVect other) {
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
  GVect& operator+=(const GVect& vect) {
    for (size_t i = 0; i < dim; ++i) {
      comp_[i] += vect.comp_[i];
    }
    return *this;
  }
  GVect& operator-=(const GVect& other) {
    for (size_t i = 0; i < dim; ++i) {
      comp_[i] -= other.comp_[i];
    }
    return *this;
  }
  GVect& operator*=(Scal k) {
    for (size_t i = 0; i < dim; ++i) {
      comp_[i] *= k;
    }
    return *this;
  }
  GVect& operator/=(Scal k) {
    for (size_t i = 0; i < dim; ++i) {
      comp_[i] /= k;
    }
    return *this;
  }
  GVect operator+(GVect other) const {
    other += *this;
    return other;
  }
  GVect operator-(GVect other) const {
    GVect tmp(*this);
    tmp -= other;
    return tmp;
  }
  GVect operator*(Scal k) const {
    GVect tmp(*this);
    tmp *= k;
    return tmp;
  }
  GVect operator/(Scal k) const {
    GVect tmp(*this);
    tmp /= k;
    return tmp;
  }
  GVect& operator*=(const GVect& other) {
    for (size_t i = 0; i < dim; ++i) {
      comp_[i] *= other.comp_[i];
    }
    return *this;
  }
  GVect operator*(const GVect& other) const {
    GVect tmp(*this);
    tmp *= other;
    return tmp;
  }
  GVect& operator/=(const GVect& other) {
    for (size_t i = 0; i < dim; ++i) {
      comp_[i] /= other.comp_[i];
    }
    return *this;
  }
  GVect operator/(const GVect& other) const {
    GVect tmp(*this);
    tmp /= other;
    return tmp;
  }
  bool operator==(const GVect& other) const {
    for (size_t i = 0; i < dim; ++i) {
      if (!(comp_[i] == other.comp_[i])) {
        return false;
      }
    }
    return true;
  }
  bool operator<(const GVect& other) const {
    for (size_t i = 0; i < dim; ++i) {
      if (!(comp_[i] < other.comp_[i])) {
        return false;
      }
    }
    return true;
  }
  bool operator<=(const GVect& other) const {
    for (size_t i = 0; i < dim; ++i) {
      if (!(comp_[i] <= other.comp_[i])) {
        return false;
      }
    }
    return true;
  }
  bool lexless(const GVect& o) const {
    return comp_ < o.comp_;
  }
  // TODO: remove, replace with GVect(0)
  static const GVect kZero;
  // TODO: remove, replace with GVect(1)
  static GVect GetUnit(size_t i) {
    GVect res = kZero;
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
  Scal dot(const GVect& other) const {
    Scal sum = 0;
    for (size_t i = 0; i < dim; ++i) {
      sum += comp_[i] * other.comp_[i];
    }
    return sum;
  }
  Scal cross_third(const GVect& other) const {
    return comp_[0] * other.comp_[1] - comp_[1] * other.comp_[0];
  }
  GVect cross(const GVect& other) const {
    const GVect& a = *this;
    const GVect& b = other;
    return GVect(a[1]*b[2]-a[2]*b[1],
                a[2]*b[0]-a[0]*b[2],
                a[0]*b[1]-a[1]*b[0]);
  }
  Scal dist(GVect other) const {
    other -= *this;
    return other.norm();
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
  // TODO: revise, may lead to undesired conversion
  template <class T=Scal>
  operator std::vector<T>() const {
    return std::vector<T>(comp_.begin(), comp_.end());
  }
  class LexLess {
   public:
    bool operator()(GVect a, GVect b) {
      return a.lexless(b);
    }
  };
};

template <class Scal, size_t dim>
const GVect<Scal, dim> GVect<Scal, dim>::kZero =
    GVect<Scal, dim>(static_cast<Scal>(0.));

template <class Scal, size_t dim>
std::ostream& operator<<(std::ostream& out, const GVect<Scal, dim>& vect) {
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
std::istream& operator>>(std::istream& in, GVect<Scal, dim>& vect) {
  for (size_t i = 0; i < dim; ++i) {
    in >> vect[i];
  }
  return in;
}

template <class _Vect>
struct Rect {
  using Vect = _Vect;
  static constexpr size_t dim = Vect::dim;

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


