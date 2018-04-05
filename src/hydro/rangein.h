#pragma once

#include "block.h"

namespace geom {

template <class Idx_, int dim_>
class GRangeIn {
  using Idx = Idx_;
  using B = GBlock<Idx, dim_>; // block 
  using I = typename B::iterator; // block iterator
  const B& ba_; // block all
  const B& bi_; // block inner

 public:
  class iterator {
    const B& ba_; // block all
    I i_;  // block iterator over bi
   public:
    explicit iterator(const B& ba, const I& i)
        : ba_(ba), i_(i)
    {}
    iterator& operator++() {
      ++i_;
      return *this;
    }
    iterator& operator--() {
      --i_;
      return *this;
    }
    bool operator==(const iterator& o /*other*/) const {
      return i_ == o.i_;
    }
    bool operator!=(const iterator& o /*other*/) const {
      return !(*this == o);
    }
    Idx operator*() const {
      return ba_.GetIdx(*i_);
    }
  };

  GRangeIn(const B& ba /*block all*/, const B& bi /*block inner*/)
      : ba_(ba), bi_(bi)
  {}
  iterator begin() const {
    return iterator(ba_, bi_.begin());
  }
  iterator end() const {
    return iterator(ba_, bi_.end());
  }
};



} // namespace geom
