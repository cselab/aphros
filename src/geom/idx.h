// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <cstddef> // ptrdiff_t

#include "vect.h"

using IntIdx = std::ptrdiff_t;

// Integer multi-index
template <size_t dim>
using GMIdx = generic::Vect<IntIdx, dim>;

// Typed index, instances distinct by id_.
template <int id_>
class GIdx {
 public:
  GIdx() : i_(0) {}

  explicit GIdx(size_t i) : i_(i) {}

  explicit operator size_t() const {
    return i_;
  }
  size_t GetRaw() const {
    return i_;
  }
  IntIdx operator-(GIdx o) const {
    return i_ - o.i_;
  }
  GIdx& operator++() {
    ++i_;
    return *this;
  }
  GIdx& operator+=(size_t add) {
    i_ += add;
    return *this;
  }
  bool operator==(GIdx o) const {
    return i_ == o.i_;
  }
  bool operator!=(GIdx o) const {
    return !(*this == o);
  }

 public:
  static constexpr int id = id_;
  size_t i_;
};

using IdxCell = GIdx<0>;
using IdxFace = GIdx<1>;
using IdxNode = GIdx<2>;
