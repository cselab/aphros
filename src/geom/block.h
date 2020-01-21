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
  using MIdx = GMIdx<dim_>;

  static constexpr size_t dim = dim_;

  class iterator {
    const GBlock* p_; // parent
    MIdx w_;

   public:
    // p: parent
    // w: current position
    iterator(const GBlock* p, MIdx w) : p_(p), w_(w) {}
    iterator& operator++() {
      for (size_t d = 0; d < dim; ++d) {
        ++w_[d];
        if (w_[d] == p_->e_[d] && d < dim - 1) {
          w_[d] = p_->b_[d];
        } else {
          break;
        }
      }
      return *this;
    }
    iterator& operator--() {
      for (size_t n = 0; n < dim; ++n) {
        if (w_[n] == p_->b_[n]) {
          w_[n] = p_->e_[n] - 1;
        } else {
          --w_[n];
          break;
        }
      }
      return *this;
    }
    bool operator==(const iterator& o) const {
      return w_ == o.w_;
    }
    bool operator!=(const iterator& o) const {
      return !(*this == o);
    }
    MIdx operator*() const {
      return w_;
    }
    // Returns the number of calls to ++ which are equivalent to ++w_[0]
    size_t GetNumLite() const {
      assert(w_[0] + 1 <= p_->GetEnd()[0]);
      return p_->GetEnd()[0] - w_[0] - 1;
    }
    // Increments iterator by `a` assuming `a <= GetLite()`
    void LiteInc(size_t a) {
      assert(a <= GetNumLite());
      w_[0] += a;
    }
  };

  GBlock(MIdx b, MIdx s) : b_(b), s_(s), e_(b_ + s_) {}
  GBlock() : GBlock(MIdx(0), MIdx(0)) {}
  GBlock(MIdx s) : GBlock(MIdx(0), s) {}
  MIdx GetBegin() const {
    return b_;
  }
  MIdx GetEnd() const {
    return e_;
  }
  MIdx GetSize() const {
    return s_;
  }
  size_t size() const {
    return s_.prod();
  }
  bool IsInside(MIdx w) const {
    return b_ <= w && w < e_;
  }
  void Clip(MIdx& w) const {
    for (size_t d = 0; d < b_.size(); ++d) {
      w[d] = std::max(b_[d], std::min(e_[d] - 1, w[d]));
    }
  }
  iterator begin() const {
    return iterator(this, b_);
  }
  iterator end() const {
    MIdx w = b_;
    w[dim - 1] = e_[dim - 1];
    return iterator(this, w);
  }

 private:
  const MIdx b_; // begin
  const MIdx s_; // size
  const MIdx e_; // end
};

// Flat index for multi-index.
template <class Idx_, size_t dim_>
class GIndex {
 public:
  using Idx = Idx_;
  using MIdx = GMIdx<dim_>;

  static constexpr size_t dim = dim_;

  // b: begin
  // s: size
  GIndex(MIdx b, MIdx s) : b_(b), s_(s), e_(b_ + s_) {}
  GIndex() : GIndex(MIdx(0), MIdx(0)) {}
  GIndex(MIdx s) : GIndex(MIdx(0), s) {}
  MIdx GetBegin() const {
    return b_;
  }
  MIdx GetEnd() const {
    return e_;
  }
  MIdx GetSize() const {
    return s_;
  }
  size_t size() const {
    return s_.prod();
  }
  GRange<Idx> Range() const {
    return GRange<Idx>(Idx(size()));
  }
  operator GRange<Idx>() const {
    return Range();
  }
  Idx GetIdx(MIdx w) const {
    w -= b_;
    size_t r = 0;
    for (size_t d = dim; d != 0;) {
      --d;
      r *= s_[d];
      r += w[d];
    }
    return Idx(r);
  }
  MIdx GetMIdx(Idx i) const {
    MIdx w;
    size_t a(i);
    for (size_t d = 0; d < dim; ++d) {
      w[d] = a % s_[d];
      a /= s_[d];
    }
    return b_ + w;
  }
  bool IsInside(MIdx w) const {
    return b_ <= w && w < e_;
  }
  void Clip(MIdx& w) const {
    for (size_t d = 0; d < b_.size(); ++d) {
      w[d] = std::max(b_[d], std::min(e_[d] - 1, w[d]));
    }
  }

 private:
  const MIdx b_; // begin
  const MIdx s_; // size
  const MIdx e_; // end
};

template <size_t dim>
using GBlockCells = GBlock<IdxCell, dim>;

template <size_t dim>
using GBlockNodes = GBlock<IdxNode, dim>;

template <size_t dim>
using GIndexCells = GIndex<IdxCell, dim>;

template <size_t dim>
using GIndexNodes = GIndex<IdxNode, dim>;
