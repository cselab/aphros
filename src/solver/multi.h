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
    if (layers.size() != size()) {
      throw std::runtime_error(
          std::string(__func__) +
          ": sizes differ, layers.size()=" + std::to_string(layers.size()) +
          " != size()=" + std::to_string(size()));
    }
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

// Returns values over stencil centered at cell c with color cl.
// Values for neighbors without color cl are filled with 0.
// sw: stencil half-width
template <class M, size_t sw>
struct GetStencil {
  static constexpr size_t sn = sw * 2 + 1;
  using Scal = typename M::Scal;

  std::array<typename M::Scal, sn * sn * sn> operator()(
      const GRange<size_t>& layers,
      const Multi<const FieldCell<typename M::Scal>*>& fc,
      const Multi<const FieldCell<typename M::Scal>*>& fccl, IdxCell c,
      typename M::Scal cl, const M& m) {
    using MIdx = typename M::MIdx;
    auto& bc = m.GetIndexCells();
    GBlock<IdxCell, M::dim> bo(MIdx(-sw), MIdx(sn));
    MIdx w = bc.GetMIdx(c);
    std::array<typename M::Scal, sn * sn * sn> uu;
    size_t k = 0;
    for (MIdx wo : bo) {
      IdxCell cn = bc.GetIdx(w + wo);
      typename M::Scal u = 0;
      for (auto j : layers) {
        if ((*fccl[j])[cn] == cl) {
          u = (*fc[j])[cn];
          break;
        }
      }
      uu[k++] = u;
    }
    return uu;
  }
};
