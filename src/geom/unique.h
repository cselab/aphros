// Created by Petr Karnakov on 21.11.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <memory>
#include <type_traits>

template <class T>
class UniquePtr {
 public:
  UniquePtr() = default;
  UniquePtr(const UniquePtr&) = delete;
  UniquePtr(UniquePtr&&) = default;

  template <class U>
  UniquePtr(UniquePtr<U>&& o) : p_(std::unique_ptr<U>(std::move(o.p_))) {}

  template <class... Args>
  UniquePtr(Args&&... args) : p_(new T(std::forward<Args>(args)...)) {}

  UniquePtr(std::nullptr_t) {}

  UniquePtr& operator=(const UniquePtr&) = delete;
  UniquePtr& operator=(UniquePtr&&) = default;

  template <class U, class... Args>
  void Set(Args&&... args) {
    p_ = std::unique_ptr<T>(new U(std::forward<Args>(args)...));
  }
  void Set(std::nullptr_t) {
    p_.reset(nullptr);
  }
  void Set(int) = delete;

  template <class U = T>
  U* Get() {
    static_assert(std::is_base_of<T, U>::value, "Not a derived class");
    return dynamic_cast<U*>(p_.get());
  }
  template <class U = T>
  const U* Get() const {
    static_assert(std::is_base_of<T, U>::value, "Not a derived class");
    return dynamic_cast<const U*>(p_.get());
  }
  T* operator->() {
    return p_.get();
  }
  const T* operator->() const {
    return p_.get();
  }

 private:
  std::unique_ptr<T> p_;
  template <class>
  friend class UniquePtr;
};
