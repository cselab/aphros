#pragma once

template <class Idx_>
class GRange {
 public:
  using Idx = Idx_;

  class iterator {
   public:
    explicit iterator(size_t i) : i_(i) {}
    iterator() = default;
    iterator(const iterator&) = default;
    iterator& operator=(const iterator&) = default;
    iterator& operator++() {
      ++i_;
      return *this;
    }
    bool operator==(const iterator& o) const {
      return i_ == o.i_;
    }
    bool operator!=(const iterator& o) const {
      return !(*this == o);
    }
    Idx operator*() const {
      return Idx(i_);
    }

   private:
    size_t i_;
  };

  GRange() : b_(0) , e_(0) {}
  GRange(size_t b, size_t e) : b_(b) , e_(e) {}

  iterator begin() const {
    return iterator(b_);
  }
  iterator end() const {
    return iterator(e_);
  }
  size_t size() const {
    return e_ - b_;
  }

 private:
  size_t b_, e_;
};


