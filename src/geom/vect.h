// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <array>
#include <cmath>
#include <istream>
#include <limits>
#include <ostream>
#include <string>
#include <vector>

template <class T>
T sqr(T a) {
  return a * a;
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

template <class Scal_, size_t dim_>
class Vect {
 public:
  static constexpr size_t dim = dim_;
  using Scal = Scal_;
  using value_type = Scal;

 public:
  friend void swap(Vect& first, Vect& second) {
    using std::swap;
    swap(first.data_, second.data_);
  }
  Vect() = default;
  Vect(const Vect& vect) = default;
  explicit Vect(Scal value) {
    for (auto& a : data_) {
      a = value;
    }
  }
  template <class... Args>
  explicit Vect(Scal first, Args... args)
      : data_{{first, static_cast<Scal>(args)...}} {
    static_assert(
        1 + sizeof...(args) == dim,
        "Vect constructor takes exactly 1 or dim arguments");
  }
  template <class T>
  explicit Vect(const Vect<T, dim>& other) {
    for (size_t i = 0; i < dim; ++i) {
      data_[i] = static_cast<Scal>(other[i]);
    }
  }
  template <class T>
  explicit Vect(const std::vector<T>& v) {
    for (size_t i = 0; i < dim; ++i) {
      data_[i] = (i < v.size() ? static_cast<Scal>(v[i]) : 0);
    }
  }
  template <class T>
  explicit Vect(const std::array<T, dim>& v) {
    for (size_t i = 0; i < dim; ++i) {
      data_[i] = static_cast<Scal>(v[i]);
    }
  }
  template <class T>
  explicit Vect(const T* v) {
    for (size_t i = 0; i < dim; ++i) {
      data_[i] = static_cast<Scal>(v[i]);
    }
  }
  template <size_t dimv>
  explicit Vect(const Vect<Scal, dimv>& v) {
    for (size_t i = 0; i < std::min(dim, dimv); ++i) {
      data_[i] = v[i];
    }
    for (size_t i = std::min(dim, dimv); i < dim; ++i) {
      data_[i] = 0;
    }
  }
  Vect& operator=(const Vect&) = default;
  constexpr size_t size() const {
    return dim;
  }
  Scal& operator[](size_t i) {
    return data_[i];
  }
  const Scal& operator[](size_t i) const {
    return data_[i];
  }
  const Scal* data() const {
    return data_.data();
  }
  Scal* data() {
    return data_.data();
  }
  Scal& back() {
    return data_.back();
  }
  const Scal& back() const {
    return data_.back();
  }
  Vect& operator+=(const Vect& vect) {
    for (size_t i = 0; i < dim; ++i) {
      data_[i] += vect.data_[i];
    }
    return *this;
  }
  Vect& operator-=(const Vect& other) {
    for (size_t i = 0; i < dim; ++i) {
      data_[i] -= other.data_[i];
    }
    return *this;
  }
  Vect& operator*=(Scal k) {
    for (size_t i = 0; i < dim; ++i) {
      data_[i] *= k;
    }
    return *this;
  }
  Vect& operator/=(Scal k) {
    for (size_t i = 0; i < dim; ++i) {
      data_[i] /= k;
    }
    return *this;
  }
  Vect& operator%=(Scal k) {
    for (size_t i = 0; i < dim; ++i) {
      data_[i] &= k;
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
  Vect operator%(Scal k) const {
    Vect tmp(*this);
    tmp %= k;
    return tmp;
  }
  Vect& operator*=(const Vect& other) {
    for (size_t i = 0; i < dim; ++i) {
      data_[i] *= other.data_[i];
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
      data_[i] /= other.data_[i];
    }
    return *this;
  }
  Vect operator/(const Vect& other) const {
    Vect tmp(*this);
    tmp /= other;
    return tmp;
  }
  Vect& operator%=(const Vect& other) {
    for (size_t i = 0; i < dim; ++i) {
      data_[i] %= other.data_[i];
    }
    return *this;
  }
  Vect operator%(const Vect& other) const {
    Vect tmp(*this);
    tmp %= other;
    return tmp;
  }
  bool operator==(const Vect& other) const {
    for (size_t i = 0; i < dim; ++i) {
      if (!(data_[i] == other.data_[i])) {
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
      if (!(data_[i] < other.data_[i])) {
        return false;
      }
    }
    return true;
  }
  bool operator>(const Vect& other) const {
    return other < (*this);
  }
  bool operator<=(const Vect& other) const {
    for (size_t i = 0; i < dim; ++i) {
      if (!(data_[i] <= other.data_[i])) {
        return false;
      }
    }
    return true;
  }
  bool operator>=(const Vect& other) const {
    return other <= (*this);
  }
  bool lexless(const Vect& o) const {
    return data_ < o.data_;
  }
  static Vect GetUnit(size_t i) {
    Vect res(0);
    res[i] = 1;
    return res;
  }
  Scal sqrnorm() const {
    Scal res = 0;
    for (size_t i = 0; i < dim; ++i) {
      res += sqr(data_[i]);
    }
    return res;
  }
  Scal norm() const {
    return std::sqrt(sqrnorm());
  }
  Scal dot(const Vect& other) const {
    Scal sum = 0;
    for (size_t i = 0; i < dim; ++i) {
      sum += data_[i] * other.data_[i];
    }
    return sum;
  }
  // Returns projection onto unit vector n.
  Vect proj(const Vect& n) const {
    return n * dot(n);
  }
  // Returns component orthogonal to unit vector n.
  Vect orth(const Vect& n) const {
    return (*this) - proj(n);
  }
  Scal cross_third(const Vect& other) const {
    return data_[0] * other.data_[1] - data_[1] * other.data_[0];
  }
  Vect cross(const Vect& other) const {
    using Vect3 = Vect<Scal, 3>;
    const Vect3 a(*this);
    const Vect3 b(other);
    return Vect(Vect3(
        a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]));
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
      r += data_[i];
    }
    return r;
  }
  Scal mean() const {
    return sum() / dim;
  }
  Scal prod() const {
    Scal r = data_[0];
    for (size_t i = 1; i < dim; ++i) {
      r *= data_[i];
    }
    return r;
  }
  Vect cumprod() const {
    Vect res(1);
    for (size_t i = 1; i < dim; ++i) {
      res[i] = res[i - 1] * data_[i - 1];
    }
    return res;
  }
  Scal norm1() const {
    Scal r = 0;
    for (size_t i = 0; i < dim; ++i) {
      r += std::abs(data_[i]);
    }
    return r;
  }
  Scal norminf() const {
    Scal r = 0;
    for (size_t i = 0; i < dim; ++i) {
      r = std::max(r, std::abs(data_[i]));
    }
    return r;
  }
  Scal max() const {
    Scal r = data_[0];
    for (size_t i = 1; i < dim; ++i) {
      r = std::max(r, data_[i]);
    }
    return r;
  }
  Scal min() const {
    Scal r = data_[0];
    for (size_t i = 1; i < dim; ++i) {
      r = std::min(r, data_[i]);
    }
    return r;
  }
  size_t argmax() const {
    size_t r = 0;
    for (size_t i = 1; i < dim; ++i) {
      if (data_[i] > data_[r]) {
        r = i;
      }
    }
    return r;
  }
  size_t argmin() const {
    size_t r = 0;
    for (size_t i = 1; i < dim; ++i) {
      if (data_[i] < data_[r]) {
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
      o[i] = std::max(data_[i], o[i]);
    }
    return o;
  }
  Vect min(Vect o) const {
    for (size_t i = 0; i < dim; ++i) {
      o[i] = std::min(data_[i], o[i]);
    }
    return o;
  }
  Vect clip(const Vect& v0, const Vect& v1) const {
    return (*this).max(v0).min(v1);
  }
  template <class T = Scal>
  operator std::vector<T>() const {
    return std::vector<T>(data_.begin(), data_.end());
  }
  template <class T = Scal>
  operator std::array<T, dim>() const {
    std::array<T, dim> r;
    for (size_t i = 0; i < dim; ++i) {
      r[i] = static_cast<T>(data_[i]);
    }
    return r;
  }
  class LexLess {
   public:
    bool operator()(Vect a, Vect b) const {
      return a.lexless(b);
    }
  };
  friend Vect operator*(Scal k, const Vect& v) {
    return v * k;
  }
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
  std::string to_string(std::string sep = " ") const {
    std::string res;
    for (size_t i = 0; i < dim; ++i) {
      if (i != 0) {
        res += sep;
      }
      res += std::to_string(data_[i]);
    }
    return res;
  }

  friend std::istream& operator>>(std::istream& in, Vect& v) {
    for (size_t i = 0; i < dim; ++i) {
      in >> v[i];
    }
    return in;
  }

 private:
  std::array<Scal, dim> data_;
};

template <class Scal, size_t dim>
const size_t Vect<Scal, dim>::dim;

} // namespace generic

template <class Scal, size_t dim>
bool IsNan(const generic::Vect<Scal, dim>& v) {
  for (size_t i = 0; i < dim; ++i) {
    if (IsNan(v[i])) {
      return true;
    }
  }
  return false;
}

template <class Scal, size_t dim>
bool IsFinite(const generic::Vect<Scal, dim>& v) {
  for (size_t i = 0; i < dim; ++i) {
    if (!IsFinite(v[i])) {
      return false;
    }
  }
  return true;
}

// Specialization providing a strict total order,
// to be used by ordered containers.
template <class T, size_t dim>
struct std::less<generic::Vect<T, dim>> {
  bool operator()(
      const generic::Vect<T, dim>& a, const generic::Vect<T, dim>& b) const {
    return a.lexless(b);
  }
};

template <class Vect>
struct Rect {
  static constexpr size_t dim = Vect::dim;
  Rect() = default;
  Rect(const Vect& low_, const Vect& high_) : low(low_), high(high_) {}
  bool IsInside(Vect x) const {
    for (size_t i = 0; i < dim; ++i) {
      if (x[i] < low[i] || high[i] < x[i]) {
        return false;
      }
    }
    return true;
  }
  Vect GetDimensions() const {
    return high - low;
  }

  Vect low, high;
};

struct GetNanHelper {
  static auto Get(double*) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  static auto Get(float*) {
    return std::numeric_limits<float>::quiet_NaN();
  }
  template <class T, size_t dim>
  static auto Get(generic::Vect<T, dim>*) {
    return generic::Vect<T, dim>(Get((T*)nullptr));
  }
};

template <class T>
T GetNan() {
  return GetNanHelper::Get((T*)nullptr);
}
