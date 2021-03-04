// Created by Petr Karnakov on 24.08.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <limits>
#include <memory>
#include <stdexcept>

#include "geom/range.h"

template <class T>
class Multi {
 public:
  Multi() = default;
  explicit Multi(size_t n) : d_(n) {}
  explicit Multi(const GRange<size_t>& layers) : d_(layers.size()) {}
  template <class... Args>
  explicit Multi(const GRange<size_t>& layers, const Args&... args)
      : d_(layers.size(), T(args...)) {}

  Multi(T u) : d_({u}) {}
  // cast to pointer
  template <class U>
  Multi(Multi<U>& u) : d_(u.size()) {
    for (size_t i = 0; i < u.size(); ++i) {
      d_[i] = &u[i];
    }
  }
  // dereference pointer
  template <class U>
  Multi(const Multi<const U*>& u) : d_(u.size()) {
    for (size_t i = 0; i < u.size(); ++i) {
      d_[i] = *u[i];
    }
  }
  // copy pointer
  template <class U>
  Multi(Multi<U*>& u) : d_(u.size()) {
    for (size_t i = 0; i < u.size(); ++i) {
      d_[i] = u[i];
    }
  }
  // cast to const
  template <class U>
  Multi(const Multi<U*>& u) : d_(u.size()) {
    for (size_t i = 0; i < u.size(); ++i) {
      d_[i] = const_cast<const U*>(u[i]);
    }
  }
  T& operator[](size_t i) {
    return d_[i];
  }
  const T& operator[](size_t i) const {
    return d_[i];
  }
  size_t size() const {
    return d_.size();
  }
  void resize(size_t n) {
    d_.resize(n);
  }
  void resize(const GRange<size_t>& layers) {
    d_.resize(layers.size());
  }
  void assert_size(const GRange<size_t>& layers) const {
    fassert_equal(layers.size(), size());
  }
  std::vector<T>& data() {
    return d_;
  }
  const std::vector<T>& data() const {
    return d_;
  }
  template <class... Args>
  void Reinit(const Args&... args) {
    for (auto& a : d_) {
      a.Reinit(args...);
    }
  }
  template <class... Args>
  void Reinit(const GRange<size_t> layers, const Args&... args) {
    resize(layers);
    Reinit(args...);
  }
  void InitAll(const T& u) {
    for (auto& a : d_) {
      a = u;
    }
  }
  std::vector<T*> GetPtr() {
    std::vector<T*> r;
    for (auto& a : d_) {
      r.push_bask(&a);
    }
    return r;
  }
  std::vector<const T*> GetConstPtr() const {
    std::vector<const T*> r;
    for (auto& a : d_) {
      r.push_bask(&a);
    }
    return r;
  }

 private:
  std::vector<T> d_;
};
