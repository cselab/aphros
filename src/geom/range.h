// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

namespace generic {

template <class T>
class Range {
 public:
  using Value = T;

  class iterator {
   public:
    explicit iterator(Value pos) : pos_(pos) {}
    iterator() = default;
    iterator(const iterator&) = default;
    iterator(iterator&&) = default;
    iterator& operator=(const iterator&) = default;
    iterator& operator=(iterator&&) = default;
    iterator& operator++() {
      ++pos_;
      return *this;
    }
    bool operator==(const iterator& other) const {
      return pos_ == other.pos_;
    }
    bool operator!=(const iterator& other) const {
      return pos_ != other.pos_;
    }
    Value operator*() const {
      return pos_;
    }

   private:
    Value pos_;
  };

  constexpr Range() {}
  constexpr explicit Range(Value end) : end_{end} {}
  constexpr Range(Value begin, Value end) : begin_{begin}, end_{end} {}

  iterator begin() const {
    return iterator(begin_);
  }
  iterator end() const {
    return iterator(end_);
  }
  constexpr size_t size() const {
    return static_cast<size_t>(end_ - begin_);
  }
  void clear() {
    (*this) = Range();
  }
  constexpr bool operator==(const Range& other) const {
    return begin_ == other.begin_ && end_ == other.end_;
  }
  constexpr bool operator!=(const Range& other) const {
    return !(*this == other);
  }

 private:
  Value begin_{0};
  Value end_{0};
};

} // namespace generic

template <class T>
using GRange = generic::Range<T>;
