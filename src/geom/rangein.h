#pragma once

#include "block.h"

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
    const B& bi_; // block inner
    I i_;  // block iterator over bi
    size_t nlite_; // number of simple inc operations
    Idx xlite_; // last idx lite
   public:
    explicit iterator(const B& ba, const B& bi, const I& i)
        : ba_(ba), bi_(bi), i_(i), nlite_(0), xlite_(ba_.GetIdx(*i_))
    {}
    iterator& operator++() {
      if (nlite_ == 0) {
        ++i_;
        nlite_ = i_.GetLite();
        xlite_ = ba_.GetIdx(*i_);
        i_.IncLite(nlite_);
      } else {
        --nlite_;
        xlite_.AddRaw(1);
      }
      return *this;
    }
    bool operator==(const iterator& o /*other*/) const {
      return xlite_ == o.xlite_;
    }
    bool operator!=(const iterator& o /*other*/) const {
      return !(*this == o);
    }
    Idx operator*() const {
      return xlite_;
    }
  };

  GRangeIn(const B& ba /*block all*/, const B& bi /*block inner*/)
      : ba_(ba), bi_(bi)
  {}
  iterator begin() const {
    return iterator(ba_, bi_, bi_.begin());
  }
  iterator end() const {
    return iterator(ba_, bi_, bi_.end());
  }
};




