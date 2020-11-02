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
