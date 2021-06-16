// Created by Petr Karnakov on 15.01.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <functional>

// analog of boost::transform_iterator
template <class Value, class Range>
class TransformIterator {
 public:
  using Iterator = decltype(((const Range*)1)->begin());
  using IterVal = decltype(**((Iterator*)1));
  using Func = std::function<Value(IterVal)>;

  class iterator {
   public:
    explicit iterator(Iterator it, const TransformIterator* owner)
        : it_(it), owner_(owner) {}
    iterator() = default;
    iterator(const iterator&) = default;
    iterator(iterator&&) = default;
    iterator& operator=(const iterator&) = default;
    iterator& operator=(iterator&&) = default;
    iterator& operator++() {
      ++it_;
      return *this;
    }
    bool operator==(const iterator& o) const {
      return it_ == o.it_;
    }
    bool operator!=(const iterator& o) const {
      return it_ != o.it_;
    }
    Value operator*() {
      return owner_->func_(*it_);
    }

   private:
    Iterator it_;
    const TransformIterator* owner_;
  };

  TransformIterator(Range range, Func func) : range_(range), func_(func) {}

  iterator begin() const {
    return iterator(range_.begin(), this);
  }
  iterator end() const {
    return iterator(range_.end(), this);
  }

 private:
  const Range range_;
  Func func_;
};

template <class V, class R>
TransformIterator<V, R> MakeTransformIterator(
    R range, typename TransformIterator<V, R>::Func func) {
  return TransformIterator<V, R>(range, func);
}
