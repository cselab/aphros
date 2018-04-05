#pragma once

#include <map>

namespace geom {

// Idx_ - instance of GIdx
template <class Value_, class Idx_>
class GMap {
 public:
  using Idx = Idx_;
  using Value = Value_;
  using Cont = std::map<size_t, Value>; // container

  GMap() {}
  explicit GMap(const GField<Value, Idx>& u) {
    for (size_t i = 0; i < u.size(); ++i) {
      d_[i] = u[Idx(i)];
    }
  }
  size_t size() const {
    return d_.size();
  }
  Value& operator[](const Idx& i) {
    return d_[i.GetRaw()];
  }
  const Value& operator[](const Idx& i) const {
    return d_[i.GetRaw()];
  }
  Value* find(const Idx& i) {
    auto it = d_.find(i.GetRaw());
    if (it != d_.end()) {
      return &it->second;
    }
    return nullptr;
  }
  const Value* find(const Idx& i) const {
    auto it = d_.find(i.GetRaw());
    if (it != d_.end()) {
      return &it->second;
    }
    return nullptr;
  }
  void erase(const Idx& i) {
    d_.erase(i.GetRaw());
  }

  class iterator {
    class Proxy {
      typename Cont::iterator it_;
     public:
      explicit Proxy(const typename Cont::iterator& it)
          : it_(it)
      {}
      Idx GetIdx() const {
        return Idx(it_->first);
      }
      Value& GetValue() {
        return it_->second;
      }
      const Value& GetValue() const {
        return it_->second;
      }
    };

    typename Cont::iterator it_;
    Proxy p_;

   public:
    explicit iterator(const typename Cont::iterator& it)
        : it_(it)
        , p_(it_)
    {}
    iterator& operator++() {
      ++it_;
      p_ = Proxy(it_);
      return *this;
    }
    iterator& operator--() {
      --it_;
      p_ = Proxy(it_);
      return *this;
    }
    bool operator==(const iterator& other) const {
      return it_ == other.it_;
    }
    bool operator!=(const iterator& other) const {
      return !(*this == other);
    }
    const Proxy& operator*() {
      return p_;
    }
    Proxy const * operator->() {
      return &p_;
    }
  };
  iterator begin() {
    return iterator(d_.begin());
  }
  iterator end() {
    return iterator(d_.end());
  }
  class const_iterator {
    class Proxy {
      typename Cont::const_iterator it_;
     public:
      explicit Proxy(const typename Cont::const_iterator& it)
          : it_(it)
      {}
      Idx GetIdx() const {
        return Idx(it_->first);
      }
      const Value& GetValue() const {
        return it_->second;
      }
    };

    typename Cont::const_iterator it_;
    Proxy p_;

   public:
    explicit const_iterator(const typename Cont::const_iterator& it)
        : it_(it)
        , p_(it_)
    {}
    const_iterator& operator++() {
      ++it_;
      p_ = Proxy(it_);
      return *this;
    }
    const_iterator& operator--() {
      --it_;
      p_ = Proxy(it_);
      return *this;
    }
    bool operator==(const const_iterator& other) const {
      return it_ == other.it_;
    }
    bool operator!=(const const_iterator& other) const {
      return !(*this == other);
    }
    const Proxy& operator*() {
      return p_;
    }
    Proxy const * operator->() {
      return &p_;
    }
  };
  const_iterator cbegin() const {
    return const_iterator(d_.cbegin());
  }
  const_iterator cend() const {
    return const_iterator(d_.cend());
  }
 
 private:
  Cont d_;
};

template <class T>
using MapCell = GMap<T, IdxCell>;

template <class T>
using MapFace = GMap<T, IdxFace>;

template <class T>
using MapNode = GMap<T, IdxNode>;


} // namespace geom
