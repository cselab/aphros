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

#if 0
// Specialization for IdxFace
// Implementation via convertion Idx->MIdx->Idx
// TODO: remove, not needed due to iterator GBlock for IdxFace
// (slow)
template <int dim>
class GRangeIn<IdxFace, dim> {
  using B = GBlockFaces<dim>; // block 
  using R = GRange<IdxFace>;  // range
  using I = typename R::iterator; // range iterator
  const B& ba_; // block all
  const B& bi_; // block inner
  R ri_; // range inner

 public:
  class iterator {
    const B& ba_; // block all
    const B& bi_; // block inner
    I i_;         // range iterator over bi
   public:
    explicit iterator(const B& ba, const B& bi, const I& i)
        : ba_(ba), bi_(bi), i_(i)
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
    IdxFace operator*() const {
      auto x = bi_.GetMIdx(*i_); 
      auto d = bi_.GetDir(*i_); 
      return ba_.GetIdx(x, d);
    }
  };

  GRangeIn(const B& ba /*block all*/, const B& bi /*block inner*/)
      : ba_(ba), bi_(bi), ri_(bi_)
  {}
  iterator begin() const {
    return iterator(ba_, bi_, ri_.begin());
  }
  iterator end() const {
    return iterator(ba_, bi_, ri_.end());
  }
};
#endif


#if 0
// Specialization for IdxFace
// Implementation with optimization: 
// don't recompute until reaching the end in x direction
// (causes segfault for Hydro)
template <int dim>
class GRangeIn<IdxFace, dim> {
  using B = GBlockFaces<dim>; // block 
  using R = GRange<IdxFace>;  // range
  using I = typename R::iterator; // range iterator
  const B& ba_; // block all
  const B& bi_; // block inner
  R ri_; // range inner

 public:
  class iterator {
    const B& ba_; // block all
    const B& bi_; // block inner
    I i_;         // range iterator over bi
    IdxFace r_;   // current result
    I in_;        // next iterator to trigger recompute of r_ via MIdx
   public:
    IdxFace IdxFromInner(I i) {
      auto x = bi_.GetMIdx(*i); 
      auto d = bi_.GetDir(*i); 
      return ba_.GetIdx(x, d); 
    }
    I GetLast(I i) {
      using Dir = typename B::Dir;
      auto x = bi_.GetMIdx(*i); 
      auto d = bi_.GetDir(*i); 
      x[0] += bi_.GetSize()[0];
      if (d != Dir::i) {
        x[0] -= 1;
      }
      return I(bi_.GetIdx(x, d).GetRaw());
    }
    explicit iterator(const B& ba, const B& bi, const I& i)
        : ba_(ba), bi_(bi), i_(i), r_(IdxFromInner(i_)), in_(i_)
    {
        in_ = GetLast(in_);
    }
    iterator& operator++() {
      if (i_ == in_) {
        ++i_;
        r_ = IdxFromInner(i_);
        ++in_;
        in_ = GetLast(in_);
      } else {
        ++i_;
        r_.AddRaw(1);
      }
      return *this;
    }
    /*
    iterator& operator--() {
      --i_;
      return *this;
    }
    */
    bool operator==(const iterator& o /*other*/) const {
      return i_ == o.i_;
    }
    bool operator!=(const iterator& o /*other*/) const {
      return !(*this == o);
    }
    IdxFace operator*() const {
      return r_;
    }
  };

  GRangeIn(const B& ba /*block all*/, const B& bi /*block inner*/)
      : ba_(ba), bi_(bi), ri_(bi_)
  {}
  iterator begin() const {
    return iterator(ba_, bi_, ri_.begin());
  }
  iterator end() const {
    return iterator(ba_, bi_, ri_.end());
  }
};
#endif


} // namespace geom
