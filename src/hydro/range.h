#pragma once



template <class _Idx>
class GRange {
  size_t pos_begin_, pos_end_;

 public:
  using Idx = _Idx;
  class iterator {
    size_t pos_;
   public:
    explicit iterator(size_t pos)
      : pos_(pos)
    {}
    iterator() = default;
    iterator(const iterator&) = default;
    iterator& operator=(const iterator&) = default;
    iterator& operator++() {
      ++pos_;
      return *this;
    }
    iterator& operator--() {
      --pos_;
      return *this;
    }
    bool operator==(const iterator& other) const {
      return pos_ == other.pos_;
    }
    bool operator!=(const iterator& other) const {
      return !(*this == other);
    }
    Idx operator*() const {
      return Idx(pos_);
    }
  };

  GRange()
    : pos_begin_(0)
    , pos_end_(0)
  {}
  GRange(size_t pos_begin, size_t pos_end)
    : pos_begin_(pos_begin)
    , pos_end_(pos_end)
  {}
  iterator begin() const {
    return iterator(pos_begin_);
  }
  iterator end() const {
    return iterator(pos_end_);
  }
  size_t size() const {
    return pos_end_ - pos_begin_;
  }
};


