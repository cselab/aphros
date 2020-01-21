// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

template <class T>
class GRange {
 public:
  using Value = T;

  class iterator {
   public:
    explicit iterator(Value i) : i_(i) {}
    iterator() = default;
    iterator(const iterator&) = default;
    iterator(iterator&&) = default;
    iterator& operator=(const iterator&) = default;
    iterator& operator=(iterator&&) = default;
    iterator& operator++() {
      ++i_;
      return *this;
    }
    bool operator==(const iterator& o) const {
      return i_ == o.i_;
    }
    bool operator!=(const iterator& o) const {
      return i_ != o.i_;
    }
    Value operator*() const {
      return i_;
    }

   private:
    Value i_;
  };

  GRange() : b_(0), e_(0) {}
  explicit GRange(Value e) : b_(0), e_(e) {}
  GRange(Value b, Value e) : b_(b), e_(e) {}

  iterator begin() const {
    return iterator(b_);
  }
  iterator end() const {
    return iterator(e_);
  }
  size_t size() const {
    return static_cast<size_t>(e_ - b_);
  }
  void clear() {
    (*this) = GRange();
  }
  bool operator==(const GRange& o) const {
    return b_ == o.b_ && e_ == o.e_;
  }
  bool operator!=(const GRange& o) const {
    return !(*this == o);
  }

 private:
  Value b_;
  Value e_;
};
