// Created by Petr Karnakov on 01.01.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <functional>

// analog of boost::filter_iterator
template <class Range>
class FilterIterator {
 public:
  using Iterator = decltype(((const Range*)1)->begin());
  using Value = decltype(**((Iterator*)1));
  using Func = std::function<bool(Value)>;

  class iterator {
   public:
    explicit iterator(Iterator it, const FilterIterator* owner)
        : it_(it), owner_(owner) {
      while (it_ != owner_->range_.end() && !owner_->func_(*it_)) {
        ++it_;
      }
    }
    iterator() = default;
    iterator(const iterator&) = default;
    iterator(iterator&&) = default;
    iterator& operator=(const iterator&) = default;
    iterator& operator=(iterator&&) = default;
    iterator& operator++() {
      do {
        ++it_;
      } while (it_ != owner_->range_.end() && !owner_->func_(*it_));
      return *this;
    }
    bool operator==(const iterator& o) const {
      return it_ == o.it_;
    }
    bool operator!=(const iterator& o) const {
      return it_ != o.it_;
    }
    Value operator*() {
      return *it_;
    }

   private:
    Iterator it_;
    const FilterIterator* owner_;
  };

  FilterIterator(Range range, Func func) : range_(range), func_(func) {}

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

template <class R>
FilterIterator<R> MakeFilterIterator(
    R range, typename FilterIterator<R>::Func func) {
  return FilterIterator<R>(range, func);
}
