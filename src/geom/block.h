// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <algorithm> // max, min
#include <cassert>

#include "idx.h"
#include "range.h"

// Block of multi-indices.
template <class Idx_, size_t dim_>
class GBlock {
 public:
  using Idx = Idx_;
  using MIdx = generic::MIdx<dim_>;

  static constexpr size_t dim = dim_;

  class iterator {
    const GBlock* owner_; // parent
    MIdx pos_;

   public:
    // p: parent
    // w: current position
    iterator(const GBlock* p, MIdx w) : owner_(p), pos_(w) {}
    iterator& operator++() {
      for (size_t d = 0; d < dim; ++d) {
        ++pos_[d];
        if (pos_[d] == owner_->end_[d] && d < dim - 1) {
          pos_[d] = owner_->begin_[d];
        } else {
          break;
        }
      }
      return *this;
    }
    iterator& operator--() {
      for (size_t n = 0; n < dim; ++n) {
        if (pos_[n] == owner_->begin_[n]) {
          pos_[n] = owner_->end_[n] - 1;
        } else {
          --pos_[n];
          break;
        }
      }
      return *this;
    }
    bool operator==(const iterator& o) const {
      return pos_ == o.pos_;
    }
    bool operator!=(const iterator& o) const {
      return !(*this == o);
    }
    MIdx operator*() const {
      return pos_;
    }
    // Returns the number of calls to ++ which are equivalent to ++pos_[0]
    size_t GetNumLite() const {
      assert(pos_[0] + 1 <= owner_->GetEnd()[0]);
      return owner_->GetEnd()[0] - pos_[0] - 1;
    }
    // Increments iterator by `a` assuming `a <= GetLite()`
    void LiteInc(size_t a) {
      assert(a <= GetNumLite());
      pos_[0] += a;
    }
  };

  GBlock(MIdx b, MIdx s) : begin_(b), size_(s), end_(begin_ + size_) {}
  GBlock() : GBlock(MIdx(0), MIdx(0)) {}
  GBlock(MIdx s) : GBlock(MIdx(0), s) {}
  MIdx GetBegin() const {
    return begin_;
  }
  MIdx GetEnd() const {
    return end_;
  }
  MIdx GetSize() const {
    return size_;
  }
  size_t size() const {
    return size_.prod();
  }
  bool IsInside(MIdx w) const {
    return begin_ <= w && w < end_;
  }
  void Clip(MIdx& w) const {
    for (size_t d = 0; d < begin_.size(); ++d) {
      w[d] = std::max(begin_[d], std::min(end_[d] - 1, w[d]));
    }
  }
  iterator begin() const {
    return iterator(this, begin_);
  }
  iterator end() const {
    MIdx w = begin_;
    w[dim - 1] = end_[dim - 1];
    return iterator(this, w);
  }

 private:
  const MIdx begin_; // begin
  const MIdx size_; // size
  const MIdx end_; // end
};

namespace generic {

template <class Idx, size_t dim>
using Block = GBlock<Idx, dim>;

} // namespace generic

// Flat index for multi-index.
template <class Idx_, size_t dim_>
class GIndex {
 public:
  using Idx = Idx_;
  using MIdx = generic::MIdx<dim_>;

  static constexpr size_t dim = dim_;

  // b: begin
  // s: size
  GIndex(MIdx b, MIdx s) : begin_(b), size_(s), end_(begin_ + size_) {}
  GIndex() : GIndex(MIdx(0), MIdx(0)) {}
  GIndex(MIdx s) : GIndex(MIdx(0), s) {}
  MIdx GetBegin() const {
    return begin_;
  }
  MIdx GetEnd() const {
    return end_;
  }
  MIdx GetSize() const {
    return size_;
  }
  size_t size() const {
    return size_.prod();
  }
  GRange<Idx> Range() const {
    return GRange<Idx>(Idx(size()));
  }
  operator GRange<Idx>() const {
    return Range();
  }
  Idx GetIdx(MIdx w) const {
    w -= begin_;
    size_t r = 0;
    for (size_t d = dim; d != 0;) {
      --d;
      r *= size_[d];
      r += w[d];
    }
    return Idx(r);
  }
  MIdx GetMIdx(Idx i) const {
    MIdx w;
    size_t a(i);
    for (size_t d = 0; d < dim; ++d) {
      w[d] = a % size_[d];
      a /= size_[d];
    }
    return begin_ + w;
  }
  bool IsInside(MIdx w) const {
    return begin_ <= w && w < end_;
  }
  void Clip(MIdx& w) const {
    for (size_t d = 0; d < begin_.size(); ++d) {
      w[d] = std::max(begin_[d], std::min(end_[d] - 1, w[d]));
    }
  }

 private:
  const MIdx begin_;
  const MIdx size_;
  const MIdx end_;
};

template <size_t dim>
using GBlockCells = GBlock<IdxCell, dim>;

template <size_t dim>
using GBlockNodes = GBlock<IdxNode, dim>;

template <size_t dim>
using GIndexCells = GIndex<IdxCell, dim>;

template <size_t dim>
using GIndexNodes = GIndex<IdxNode, dim>;
