#pragma once

#include <functional>

template <class Range>
class Filter {
 public:
  using Iterator = decltype(((const Range*)0)->begin());
  using Value = decltype(*Iterator());
  using Func = std::function<bool(Value)>;

  class iterator {
   public:
    explicit iterator(Iterator it, const Filter* owner)
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
    const Filter* owner_;
  };

  Filter(Range range, Func func) : range_(range), func_(func) {}

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
Filter<R> MakeFilter(R range, typename Filter<R>::Func func) {
  return Filter<R>(range, func);
}
