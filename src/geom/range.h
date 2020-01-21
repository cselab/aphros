// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

template <class Iter_>
class GRange {
 public:
  using Iter = Iter_;

  class iterator {
   public:
    explicit iterator(Iter it) : it_(it) {}
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
    Iter operator*() const {
      return it_;
    }

   private:
    Iter it_;
  };

  GRange() : b_(0), e_(0) {}
  explicit GRange(Iter e) : b_(0), e_(e) {}
  GRange(Iter b, Iter e) : b_(b), e_(e) {}

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
  Iter b_;
  Iter e_;
};
