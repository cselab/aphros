// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <cstddef> // ptrdiff_t

#include "vect.h"

using IntIdx = std::ptrdiff_t;

namespace generic {

template <size_t dim>
using MIdx = generic::Vect<IntIdx, dim>;

template <int id_>
class Idx {
 public:
  Idx() : value_(0) {}
  explicit Idx(size_t value) : value_(value) {}
  explicit operator size_t() const {
    return value_;
  }
  size_t GetRaw() const {
    return value_;
  }
  size_t raw() const {
    return value_;
  }
  IntIdx operator-(Idx o) const {
    return value_ - o.value_;
  }
  Idx& operator++() {
    ++value_;
    return *this;
  }
  Idx& operator+=(size_t add) {
    value_ += add;
    return *this;
  }
  Idx operator+(size_t add) const {
    return Idx(value_ + add);
  }
  bool operator==(Idx o) const {
    return value_ == o.value_;
  }
  bool operator!=(Idx o) const {
    return !(*this == o);
  }
  bool operator<(Idx o) const {
    return value_ < o.value_;
  }

 public:
  static constexpr int id = id_;
  size_t value_;
};

} // namespace generic

using IdxCell = generic::Idx<0>;
using IdxFace = generic::Idx<1>;
using IdxNode = generic::Idx<2>;
using IdxNci = generic::Idx<3>; // neighbor cell or face index
