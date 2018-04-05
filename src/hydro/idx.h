#pragma once

#include <cstddef> // ptrdiff_t

#include "vect.hpp"

namespace geom {

using IntIdx = std::ptrdiff_t;

// Integer multi-index
template <size_t dim>
using GMIdx = GVect<IntIdx, dim>;

// Typed integer index, instances distinct by id_.
template <int id_>
class GIdx {
  static constexpr int id = id_; 
  size_t i_;
  static constexpr size_t kNone = -1;
 public:
  GIdx() {}
  explicit GIdx(size_t i)
    : i_(i)
  {}
  size_t GetRaw() const {
    return i_;
  }
  void AddRaw(IntIdx add) {
    i_ += add;
  }
  bool operator==(GIdx o) const {
    return i_ == o.i_;
  }
  bool operator!=(GIdx o) const {
    return !(*this == o);
  }
  static GIdx None() {
    return GIdx(-1);
  }
  bool IsNone() const {
    return *this == None();
  }
};

} // namespace geom
