// Created by Petr Karnakov on 23.07.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <utility>

#include "block.h"
#include "dir.h"

// Specialization for IdxFace, d-major
// Enumeration order (to more frequent): z, y, x, dir
template <size_t dim_>
class GBlock<IdxFace, dim_> { // [s]andbox
 public:
  static constexpr size_t dim = dim_;
  using Idx = IdxFace;
  using MIdx = generic::MIdx<dim>;
  using Dir = GDir<dim>;

  static_assert(dim == 3, "GBlock<IdxFace,...> implemented only for dim=3");

  class iterator {
    const GBlock* o_; // owner
    MIdx x_;
    Dir d_;

   public:
    explicit iterator(const GBlock* o, MIdx x, Dir d) : o_(o), x_(x), d_(d) {}
    iterator& operator++() {
      auto& s = o_->cs_;
      auto& x = x_;
      auto& b = o_->b_;
      auto& e = o_->e_;
      auto& d = d_;

      if (x[2] == e[2]) { // z-boundary
        ++x[0];
        if (x[0] == e[0]) {
          x[0] = b[0];
          ++x[1]; // end would be: (b[0], e[1], e[2])
        }
      } else if (x[1] == e[1]) { // y-edge
        ++x[0];
        if (x[0] == e[0]) { // next x-end
          ++x[2];
          x[1] = b[1];
          x[0] = b[0];
          d = Dir(0);
          if (x[2] == e[2]) { // next z-edge
            d = Dir(2);
          }
        }
      } else if (x[0] == e[0]) { // x-edge
        x[0] = b[0];
        ++x[1];
        if (x[1] == e[1]) { // next y-edge
          d = Dir(1);
        }
      } else if (d == Dir(2)) { // dirz
        ++x[0];
        d = Dir(0);
      } else {
        d = Dir(size_t(d) + 1);
      }

      return *this;
    }
    bool operator==(const iterator& o) const {
      return x_ == o.x_ && d_ == o.d_;
    }
    bool operator!=(const iterator& o) const {
      return !(*this == o);
    }
    std::pair<MIdx, Dir> operator*() const {
      return std::make_pair(x_, d_);
    }
  };

  GBlock(MIdx begin, MIdx cells) : b_(begin), cs_(cells), e_(b_ + cs_) {}
  GBlock(MIdx cells) : GBlock(MIdx(0), cells) {}
  GBlock() : GBlock(MIdx(0), MIdx(0)) {}
  size_t size() const {
    size_t res = 0;
    for (size_t i = 0; i < dim; ++i) {
      res += GetNumFaces(i);
    }
    return res;
  }
  operator GRange<Idx>() const {
    return GRange<Idx>(0, size());
  }
  Idx GetIdx(MIdx x, Dir d) const {
    x -= b_;
    size_t r = 0;
    auto& s = cs_;
    r += (s[1] * s[0] + (s[1] + 1) * s[0] + (s[0] + 1) * s[1]) * x[2];
    if (x[2] == s[2]) {
      r += s[0] * x[1] + x[0]; // z-boundary
    } else {
      r += (3 * s[0] + 1) * x[1];
      if (x[1] == s[1]) {
        r += x[0]; // y-boundary
      } else {
        r += 3 * x[0] + size_t(d);
      }
    }
    return Idx(r);
  }
  Idx GetIdx(std::pair<MIdx, Dir> p) const {
    return GetIdx(p.first, p.second);
  }
  MIdx GetMIdx(Idx idx) const {
    return GetMIdxDir(idx).first;
  }
  Dir GetDir(Idx idx) const {
    return GetMIdxDir(idx).second;
  }
  MIdx GetSize() const {
    return cs_;
  }
  iterator begin() const {
    return iterator(this, b_, Dir(0));
  }
  iterator end() const {
    return iterator(this, MIdx(b_[0], e_[1], e_[2]), Dir(2));
  }
  std::pair<MIdx, Dir> GetMIdxDir(Idx i) const {
    size_t r = i.GetRaw();
    MIdx x;
    Dir d;
    auto& s = cs_;
    const size_t n = size();
    if (r >= n - s[0] * s[1]) { // z-boundary
      r -= n - s[0] * s[1];
      x[2] = s[2];
      x[1] = r / s[0];
      x[0] = r % s[0];
      d = Dir(2);
    } else {
      auto w = (s[1] * s[0] + (s[1] + 1) * s[0] + (s[0] + 1) * s[1]);
      x[2] = r / w;
      r = r % w;
      if (r >= w - s[0]) { // y-boundary
        r -= w - s[0];
        x[1] = s[1];
        x[0] = r;
        d = Dir(1);
      } else {
        auto q = 3 * s[0] + 1;
        x[1] = r / q;
        r = r % q;
        x[0] = r / 3;
        d = Dir(r % 3);
      }
    }
    return std::make_pair(b_ + x, d);
  }

 private:
  MIdx b_; // begin
  MIdx cs_; // cells size
  MIdx e_; // end
  size_t GetNumFaces(Dir dir) const {
    MIdx bcs = cs_;
    ++bcs[size_t(dir)];
    size_t res = 1;
    for (size_t i = 0; i < dim; ++i) {
      res *= bcs[i];
    }
    return res;
  }
  size_t GetNumFaces(size_t i) const {
    return GetNumFaces(Dir(i));
  }
};
