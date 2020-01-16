// Created by Petr Karnakov on 30.05.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <algorithm> // min, max
#include <cmath> // abs
#include <iostream> // cerr, scientific, setprecision

template <class T>
bool Cmp(T a, T b) {
  return a == b;
}

template <>
bool Cmp<double>(double a, double b) {
  return std::abs(a - b) < 1e-10;
}

template <class Idx, class B, class Scal>
Scal DiffMax(
    const B& b, const GField<Scal, Idx>& u, const GField<Scal, Idx>& v,
    const GField<bool, Idx>& mask) {
  Scal r = 0;
  for (auto i : GRange<Idx>(b)) {
    if (mask[i]) {
      r = std::max(r, std::abs(u[i] - v[i]));
    }
  }
  return r;
}

template <class Idx, class B, class Scal>
Scal Max(
    const B& b, const GField<Scal, Idx>& u, const GField<bool, Idx>& mask) {
  Scal r = 0;
  for (auto i : GRange<Idx>(b)) {
    if (mask[i]) {
      r = std::max(r, u[i]);
    }
  }
  return r;
}

template <class Idx, class B, class Scal>
Scal Mean(
    const B& b, const GField<Scal, Idx>& u, const GField<bool, Idx>& mask) {
  Scal r = 0;
  Scal w = 0.;
  for (auto i : GRange<Idx>(b)) {
    if (mask[i]) {
      r += u[i];
      w += 1.;
    }
  }
  return r / w;
}

template <class Idx, class M>
typename M::Scal DiffMax(
    const GField<typename M::Scal, Idx>& u,
    const GField<typename M::Scal, Idx>& v, const M& m,
    const GField<bool, Idx>& mask) {
  using Scal = typename M::Scal;
  Scal r = 0;
  for (auto i : m.template Get<Idx>()) {
    if (mask[i]) {
      r = std::max(r, std::abs(u[i] - v[i]));
    }
  }
  return r;
}

template <class Idx, class M>
typename M::Scal Max(
    const GField<typename M::Scal, Idx>& u, const M& m,
    const GField<bool, Idx>& mask) {
  using Scal = typename M::Scal;
  Scal r = 0;
  for (auto i : m.template Get<Idx>()) {
    if (mask[i]) {
      r = std::max(r, u[i]);
    }
  }
  return r;
}

template <class Idx, class M>
typename M::Scal Mean(
    const GField<typename M::Scal, Idx>& u, const M& m,
    const GField<bool, Idx>& mask) {
  using Scal = typename M::Scal;
  Scal r = 0;
  Scal w = 0.;
  for (auto i : m.template Get<Idx>()) {
    if (mask[i]) {
      r += u[i];
      w += 1.;
    }
  }
  return r / w;
}

#define CMP(a, b) assert(Cmp(a, b));

// Print CMP
#define PCMP(a, b)                                                        \
  std::cerr << std::scientific << std::setprecision(16) << #a << "=" << a \
            << ", " << #b << "=" << b << std::endl;                       \
  CMP(a, b);

// Print CMP if false
#define PFCMP(a, b, ftl)                                                \
  if (!Cmp(a, b)) {                                                     \
    std::cerr << std::scientific << std::setprecision(16)               \
              << "Failed cmp: " << std::endl                            \
              << #a << "=" << a << ", " << #b << "=" << b << std::endl; \
    assert(!ftl);                                                       \
  }
