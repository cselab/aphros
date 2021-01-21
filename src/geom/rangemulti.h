// Created by Petr Karnakov on 02.11.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <array>

namespace generic {

template <class Idx_, size_t dim_, size_t step_>
class RangeMulti {
 public:
  static constexpr size_t dim = dim_;
  using Idx = Idx_;
  using MIdx = std::array<int, dim>;
  static constexpr size_t step = step_;

  class iterator {
   public:
    // pos: current index between MIdx(0) and size-1
    // flat_begin: flat index of corresponding to pos=MIdx(0)
    // size: size of subblock
    // lead: size of whole block
    explicit iterator(
        const MIdx& pos, Idx flat_begin, const MIdx& size, const MIdx& lead)
        : size_(size), pos_(pos) {
      flat_ = Idx(flat_begin);
      size_t prod = 1;
      for (size_t i = 0; i < dim; ++i) {
        flat_ += prod * pos[i];
        prodone_[i] = prod;
        prod *= lead[i];
        prod_[i] = prod;
      }
    }
    iterator& operator++() {
      flat_ += step;
      for (size_t d = 0; d < dim; ++d) {
        pos_[d] += (d == 0 ? step : 1);
        if (pos_[d] >= size_[d] && d + 1 < dim) {
          flat_ += prod_[d] - prodone_[d] * pos_[d];
          pos_[d] = 0;
        } else {
          break;
        }
      }
      return *this;
    }
    bool operator==(const iterator& other) const {
      return flat_ == other.flat_;
    }
    bool operator!=(const iterator& other) const {
      return !(*this == other);
    }
    Idx operator*() const {
      return flat_;
    }

   private:
    const MIdx size_;
    MIdx offset_;
    MIdx prodone_;
    MIdx prod_;
    MIdx pos_;
    Idx flat_;
  };

  RangeMulti(Idx flat_begin, const MIdx& size, const MIdx& lead)
      : flat_begin_(flat_begin), size_(size), lead_(lead) {
    static_assert(step >= 1, "");
    assert(size <= lead);
  }
  iterator begin() const {
    MIdx first;
    for (auto& a : first) {
      a = 0;
    }
    return iterator(first, flat_begin_, size_, lead_);
  }
  iterator end() const {
    MIdx last;
    for (auto& a : last) {
      a = 0;
    }
    last[dim - 1] = size_[dim - 1];
    return iterator(last, flat_begin_, size_, lead_);
  }

 private:
  const Idx flat_begin_;
  const MIdx size_;
  const MIdx lead_;
};

} // namespace generic
