// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include "block.h"

template <class Idx_, int dim_>
class GRangeIn {
 public:
  using Idx = Idx_;
  using Block = GBlock<Idx, dim_>;
  using Indexer = GIndex<Idx, dim_>;

  class iterator {
   public:
    explicit iterator(
        const Indexer& indexer, const Block& block,
        const typename Block::iterator& block_iter)
        : indexer_(indexer)
        , block_(block)
        , block_iter_(block_iter)
        , nlite_(0)
        , idx_(indexer_.GetIdx(*block_iter_)) {}
    iterator& operator++() {
      if (nlite_ == 0) {
        ++block_iter_;
        nlite_ = block_iter_.GetNumLite();
        idx_ = indexer_.GetIdx(*block_iter_);
        block_iter_.LiteInc(nlite_);
      } else {
        --nlite_;
        ++idx_;
      }
      return *this;
    }
    bool operator==(const iterator& o) const {
      return idx_ == o.idx_;
    }
    bool operator!=(const iterator& o) const {
      return !(*this == o);
    }
    Idx operator*() const {
      return idx_;
    }

   private:
    const Indexer& indexer_;
    const Block& block_;
    typename Block::iterator block_iter_;
    size_t nlite_; // number of calls of operator++ equivalent to ++idx_
    Idx idx_;
  };

  GRangeIn(const Indexer& indexer, const Block& block)
      : indexer_(indexer), block_(block) {}
  iterator begin() const {
    return iterator(indexer_, block_, block_.begin());
  }
  iterator end() const {
    return iterator(indexer_, block_, block_.end());
  }

 private:
  const Indexer& indexer_;
  const Block& block_;
};

template <class Idx_, size_t dim_, size_t step_>
class GRangeMulti {
  using Idx = Idx_;
  using MIdx = GMIdx<dim_>;
  static constexpr size_t dim = dim_;
  static constexpr size_t step = step_;

 public:
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

  GRangeMulti(Idx flat_begin, const MIdx& size, const MIdx& lead)
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
