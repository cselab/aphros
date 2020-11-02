// Created by Petr Karnakov on 02.11.2020
// Copyright 2020 ETH Zurich

#pragma once

#include "idx.h"

namespace generic {

template <class Idx_, size_t dim_, size_t step_>
class RangeMulti {
 public:
  using Idx = Idx_;
  using MIdx = GMIdx<dim_>;
  static constexpr size_t dim = dim_;
  static constexpr size_t step = step_;

  class iterator {
   public:
    // pos: current index between MIdx(0) and size-1
    // flat_begin: flat index of corresponding to pos=MIdx(0)
    // size: size of subblock
    // lead: size of whole block
    explicit iterator(
        const MIdx& pos, Idx flat_begin, const MIdx& size, const MIdx& lead)
        : size_(size)
        , offset_(Product(lead) - size * ProductOne(lead))
        , pos_(pos)
        , flat_(size_t(flat_begin) + (pos * ProductOne(lead)).sum()) {}
    iterator& operator++() {
      flat_ += step;
      for (size_t d = 0; d < dim; ++d) {
        pos_[d] += (d == 0 ? step : 1);
        if (pos_[d] == size_[d] && d + 1 < dim) {
          pos_[d] = 0;
          flat_ += offset_[d];
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
    // 1, size[0], size[0] * size[1], ...
    static MIdx ProductOne(MIdx size) {
      MIdx res(1);
      for (size_t d = 1; d < dim; ++d) {
        res[d] = res[d - 1] * size[d - 1];
      }
      return res;
    }
    // size[0], size[0] * size[1], size[0] * size[1] * size[2]
    static MIdx Product(MIdx size) {
      MIdx res(size);
      for (size_t d = 1; d < dim; ++d) {
        res[d] = res[d - 1] * size[d];
      }
      return res;
    }

    const MIdx size_;
    const MIdx offset_; // offset in flat index after reaching size
    MIdx pos_;
    Idx flat_;
  };

  RangeMulti(Idx flat_begin, const MIdx& size, const MIdx& lead)
      : flat_begin_(flat_begin), size_(size), lead_(lead) {
    assert(size <= lead);
    assert(step >= 1);
    assert(size[0] % step == 0);
  }
  iterator begin() const {
    return iterator(MIdx(0), flat_begin_, size_, lead_);
  }
  iterator end() const {
    MIdx last(0);
    last[dim - 1] = size_[dim - 1];
    return iterator(last, flat_begin_, size_, lead_);
  }

 private:
  const Idx flat_begin_;
  const MIdx size_;
  const MIdx lead_;
};

} // namespace generic
