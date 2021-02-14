// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <utility>

#include "block.h"
#include "dir.h"

// Specialization for IdxFace
// Enumeration order (to more frequent): dir, z, y, x
template <size_t dim_>
class GBlock<IdxFace, dim_> {
 public:
  static constexpr size_t dim = dim_;
  using Idx = IdxFace;
  using MIdx = generic::MIdx<dim>;
  using Dir = GDir<dim>;

  class iterator {
    const GBlock* owner_;
    MIdx w_;
    Dir dir_;

   public:
    explicit iterator(const GBlock* owner, MIdx w, Dir d)
        : owner_(owner), w_(w), dir_(d) {}
    iterator& operator++() {
      MIdx wd(dir_);
      for (size_t n = 0; n < dim; ++n) {
        ++w_[n]; // increment current
        if (w_[n] ==
            owner_->start_[n] + owner_->ncells_[n] + wd[n]) { // if end reached
          if (n < dim - 1) {
            w_[n] = owner_->start_[n]; // reset to begin
          } else { // if end reached for last dim
            if (size_t(dir_) < dim - 1) {
              dir_ = Dir(size_t(dir_) + 1);
              w_ = owner_->start_;
            }
          }
        } else {
          break;
        }
      }
      return *this;
    }
    bool operator==(const iterator& o) const {
      return w_ == o.w_ && dir_ == o.dir_;
    }
    bool operator!=(const iterator& o) const {
      return !(*this == o);
    }
    std::pair<MIdx, Dir> operator*() const {
      return std::make_pair(w_, dir_);
    }
    // Returns number of calls ++i_ for which
    // GetIdx(*++i_).GetRaw() == GetIdx(*i_).GetRaw() + 1
    // (i.e. increment of iterator equivalent to increment of index)
    size_t GetNumLite() const {
      return owner_->start_[0] + owner_->ncells_[0] +
             (dir_.raw() == 0 ? 1 : 0) - w_[0] - 1;
    }
    void LiteInc(size_t a) {
      assert(a <= GetNumLite());
      w_[0] += a;
    }
  };

  // start: begin, lower corner cell index
  // ncells: cells size
  GBlock(MIdx start, MIdx ncells)
      : start_(start)
      , ncells_(ncells)
      , count_(MIdx(ncells_.prod()) + MIdx(ncells_.prod()) / ncells_) {}
  GBlock(MIdx ncells) : GBlock(MIdx(0), ncells) {}
  GBlock() : GBlock(MIdx(0), MIdx(0)) {}
  size_t size() const {
    return count_.sum();
  }
  operator GRange<Idx>() const {
    return GRange<Idx>(0, size());
  }
  MIdx GetSize() const {
    return ncells_;
  }
  iterator begin() const {
    return iterator(this, start_, Dir(0));
  }
  iterator end() const {
    MIdx w = start_;
    w[dim - 1] = (start_ + ncells_)[dim - 1] + 1;
    return iterator(this, w, Dir(dim - 1));
  }

 private:
  const MIdx start_;
  const MIdx ncells_; // cells size
  const MIdx count_; // number of faces in each direction
};

// Specialization for IdxFace
// Enumeration order (to more frequent): dir, z, y, x
template <size_t dim_>
class GIndex<IdxFace, dim_> {
 public:
  static constexpr size_t dim = dim_;
  using Idx = IdxFace;
  using MIdx = generic::MIdx<dim>;
  using Dir = GDir<dim>;

  // start: begin, lower corner
  // rawsize: raw size (underlying block that fits all faces_)
  GIndex(MIdx start, MIdx rawsize)
      : indexc_(start, rawsize), ncells_(indexc_.size()) {}
  GIndex(MIdx rawsize) : GIndex(MIdx(0), rawsize) {}
  GIndex() : GIndex(MIdx(0), MIdx(0)) {}
  size_t size() const {
    return ncells_ * dim;
  }
  operator GRange<Idx>() const {
    return GRange<Idx>(0, size());
  }
  Idx GetIdx(MIdx w, Dir d) const {
    return Idx(d.raw() * ncells_ + indexc_.GetIdx(w));
  }
  Idx GetIdx(std::pair<MIdx, Dir> p) const {
    return GetIdx(p.first, p.second);
  }
  Dir GetDir(Idx f) const {
    const auto r = f.raw();
    const auto n = ncells_;
    for (size_t d = 0; d < dim; ++d) {
      if (r < n * (d + 1)) {
        return Dir(d);
      }
    }
    return {};
  }
  std::pair<MIdx, Dir> GetMIdxDir(Idx f) const {
    const auto r = f.raw();
    const auto n = ncells_;
    for (size_t d = 0; d < dim; ++d) {
      if (r < n * (d + 1)) {
        return {indexc_.GetMIdx(r - n * d), Dir(d)};
      }
    }
    return {};
  }
  MIdx GetMIdx(Idx f) const {
    const auto r = f.raw();
    const auto n = ncells_;
    for (size_t d = 0; d < dim; ++d) {
      if (r < n * (d + 1)) {
        return indexc_.GetMIdx(r - n * d);
      }
    }
    return {};
  }

 private:
  const GIndex<size_t, dim> indexc_; // indexer over cells
  const size_t ncells_; // total number of cells with padding
};

template <size_t dim>
using GBlockFaces = GBlock<IdxFace, dim>;

template <size_t dim>
using GIndexFaces = GIndex<IdxFace, dim>;
