// Created by Petr Karnakov on 30.10.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <vector>

namespace fast_allocator_impl {

class FixedAllocatorGeneric {
 public:
  virtual ~FixedAllocatorGeneric() {}
  virtual void* Allocate() = 0;
  virtual void Deallocate(void* ptr) noexcept = 0;
};

template <size_t ChunkSize, size_t PoolSize>
class FixedAllocator : public FixedAllocatorGeneric {
  using pointer = char*;
  using size_type = std::size_t;
  static constexpr size_type kEffChunkSize =
      (ChunkSize > sizeof(pointer) ? ChunkSize : sizeof(pointer));

  pointer head_;
  std::vector<pointer> vpool_;

  pointer AllocatePool();

 public:
  FixedAllocator() : head_(nullptr) {}
  ~FixedAllocator();
  void* Allocate() override;
  void Deallocate(void* ptr) noexcept override;
};

template <size_t ChunkSize, size_t PoolSize>
auto FixedAllocator<ChunkSize, PoolSize>::AllocatePool() -> pointer {
  vpool_.reserve(vpool_.size() + 1);
  pointer begin = static_cast<pointer>(::operator new(PoolSize));
  vpool_.push_back(begin);
  pointer end = begin + PoolSize;
  pointer pos;
  for (pos = begin; pos + kEffChunkSize < end; pos += kEffChunkSize) {
    *static_cast<pointer*>(static_cast<void*>(pos)) = pos + kEffChunkSize;
  }
  if (pos < end) {
    *static_cast<pointer*>(static_cast<void*>(pos)) = nullptr;
  }
  return begin;
}

template <size_t ChunkSize, size_t PoolSize>
void* FixedAllocator<ChunkSize, PoolSize>::Allocate() {
  if (!head_) {
    head_ = AllocatePool();
  }

  pointer res = head_;
  head_ = *static_cast<pointer*>(static_cast<void*>(head_));
  return res;
}

template <size_t ChunkSize, size_t PoolSize>
void FixedAllocator<ChunkSize, PoolSize>::Deallocate(void* ptr) noexcept {
  *static_cast<pointer*>(ptr) = head_;
  head_ = static_cast<pointer>(ptr);
}

template <size_t ChunkSize, size_t PoolSize>
FixedAllocator<ChunkSize, PoolSize>::~FixedAllocator() {
  for (auto pool_ptr : vpool_) {
    ::operator delete(pool_ptr);
  }
}

extern std::map<size_t, FixedAllocatorGeneric*> g_alloc_map;

// Contains a series of FixedAllocator objects.
// The chunk size starting from ChunkSize and being doubled on each level.
template <size_t ChunkSize, size_t Depth, size_t PoolSize = (1 << 16)>
class FixedAllocatorSeries {
  FixedAllocatorSeries<ChunkSize * 2, Depth - 1, PoolSize> child_;
  FixedAllocator<ChunkSize, PoolSize> alloc_;

 public:
  FixedAllocatorSeries() {
    g_alloc_map[ChunkSize] = &alloc_;
  }
};

template <size_t ChunkSize, size_t PoolSize>
class FixedAllocatorSeries<ChunkSize, 0, PoolSize> {
 public:
  FixedAllocatorSeries() {}
};

template <class T>
class FastAllocator {
 public:
  using value_type = T;
  using pointer = T*;
  using const_pointer = const T*;
  using reference = T&;
  using const_reference = const T&;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  template <class U>
  struct rebind {
    using other = FastAllocator<U>;
  };
  pointer allocate(size_type n);
  void deallocate(pointer p, size_type n);
  void construct(pointer p, const_reference val);
  void destroy(pointer p);
  size_type max_size() const;
  pointer address(reference x) const;
  const_pointer address(const_reference x) const;

  bool operator==(const FastAllocator& other) {
    return true;
  }
  bool operator!=(const FastAllocator& other) {
    return false;
  }

  FastAllocator() {}

 private:
  size_type GetNearestRegisteredSize(size_type n) {
    auto pos = g_alloc_map.lower_bound(n);
    return pos == g_alloc_map.end() ? 0 : pos->first;
  }
};

template <class T>
auto FastAllocator<T>::allocate(size_type n) -> pointer {
  size_type bytes = n * sizeof(T);
  size_type nearest_bytes = GetNearestRegisteredSize(bytes);
  if (nearest_bytes) {
    return static_cast<pointer>(g_alloc_map[nearest_bytes]->Allocate());
  }
  return static_cast<pointer>(operator new(bytes));
}

template <class T>
void FastAllocator<T>::deallocate(pointer p, size_type n) {
  size_type bytes = n * sizeof(T);
  size_type nearest_bytes = GetNearestRegisteredSize(bytes);
  if (nearest_bytes) {
    return g_alloc_map[nearest_bytes]->Deallocate(p);
  }
  return operator delete(p);
}

template <class T>
void FastAllocator<T>::construct(pointer p, const_reference val) {
  new (static_cast<void*>(p)) T(val);
}

template <class T>
void FastAllocator<T>::destroy(pointer p) {
  p->~T();
}

template <class T>
auto FastAllocator<T>::max_size() const -> size_type {
  return std::numeric_limits<size_type>::max();
}

template <class T>
auto FastAllocator<T>::address(reference x) const -> pointer {
  return std::addressof(x);
}

template <class T>
auto FastAllocator<T>::address(const_reference x) const -> const_pointer {
  return std::addressof(x);
}

} // namespace fast_allocator_impl

template <size_t ChunkSize, size_t PoolSize>
using FixedAllocator = fast_allocator_impl::FixedAllocator<ChunkSize, PoolSize>;

template <class T>
using FastAllocator = fast_allocator_impl::FastAllocator<T>;

template <size_t ChunkSizeInitial, size_t Depth, size_t PoolSize>
using FixedAllocatorSeries = fast_allocator_impl::FixedAllocatorSeries<
    ChunkSizeInitial, Depth, PoolSize>;
