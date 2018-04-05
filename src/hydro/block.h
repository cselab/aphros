#pragma once

#include <algorithm> // max, min

#include "idx.h"
#include "range.h"

namespace geom {

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
    explicit iterator(const GBlock* p /*parent*/, MIdx w)
      : p_(p), w_(w)
    {}
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
  };

  GBlock()
    : b_(0), s_(0), e_(0)
  {}
  GBlock(MIdx s /*size*/)
    : b_(0), s_(s), e_(b_ + s_)
  {}
  GBlock(MIdx b /*begin*/, MIdx s /*size*/)
    : b_(b), s_(s), e_(b_ + s_)
  {}
  // TODO: rename to GetSize
  MIdx GetDimensions() const {
    return s_;
  }
  MIdx GetBegin() const {
    return b_;
  }
  MIdx GetEnd() const {
    return e_;
  }
  size_t size() const {
    size_t r = 1;
    for (size_t d = 0; d < dim; ++d) {
      r *= s_[d];
    }
    return r;
  }
  operator GRange<Idx>() const {
    return GRange<Idx>(0, size());
  }
  Idx GetIdx(MIdx w) const {
    w -= b_;
    size_t r = 0;
    for (size_t d = dim; d != 0; ) {
      --d;
      r *= s_[d];
      r += w[d];
    }
    return Idx(r);
  }
  MIdx GetMIdx(Idx i) const {
    MIdx w;
    size_t a = i.GetRaw();
    for (size_t d = 0; d < dim; ++d) {
      w[d] = a % s_[d];
      a /= s_[d];
    }
    return b_ + w;
  }
  bool IsInside(MIdx w) const {
    return b_ <= w && w < e_;
  }
  // TODO: rename to clip
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
  MIdx b_; // begin
  MIdx s_; // size
  MIdx e_; // end
};

template <size_t dim>
using GBlockCells = GBlock<IdxCell, dim>;

template <size_t dim>
using GBlockNodes = GBlock<IdxNode, dim>;

} // namespace geom
