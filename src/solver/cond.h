// Created by Petr Karnakov on 15.05.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <memory>
#include <string>
#include <utility>

#include "geom/map.h"

class CondCell {
 public:
  virtual ~CondCell() {}
};

// Condition for value
template <class V>
class CondCellVal : public CondCell {
 public:
  virtual V second() const = 0;
};

// Given value
template <class V>
class CondCellValFixed : public CondCellVal<V> {
 public:
  explicit CondCellValFixed(const V& v) : v_(v) {}
  V second() const override {
    return v_;
  }
  void Set(const V& v) {
    v_ = v;
  }

 private:
  V v_;
};

enum class BCondType {
  dirichlet, // u = val
  neumann, // du/dn = val
  mixed, // if u is vector,
         //   u.proj(n) = val.proj(n)
         //   (du/dn).orth(n) = val.orth(n)
         // if u is scalar, equivalent to `dirichlet`
  reflect, // equivalent to `mixed` with val=0
  extrap, // extrapolated from neighbor cells
};

// Boundary condition on face or embedded boundary.
template <class T>
struct BCond {
  BCond() = default;
  BCond(BCondType type_, size_t nci_, T val_)
      : type(type_), nci(nci_), val(val_) {}
  BCond(BCondType type_, size_t nci_) : type(type_), nci(nci_) {}
  BCondType type = BCondType::dirichlet;
  size_t nci = 0; // neighbor cell id on faces and 0 on embedded boundaries
  T val = T(0); // field value (dirichlet) or normal gradient (neumann)
};
