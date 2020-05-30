// Created by Petr Karnakov on 15.05.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <memory>
#include <string>
#include <utility>

#include "geom/map.h"
#include "geom/unique.h"

class CondFace {
 public:
  virtual ~CondFace() {}
  virtual size_t GetNci() const {
    return nci_;
  }

 protected:
  // nci: neighbour cell id
  CondFace(size_t nci) : nci_(nci) {}

 private:
  size_t nci_;
};

// First order exptrapolation from cell
class CondFaceExtrap : public CondFace {
 public:
  CondFaceExtrap(size_t nci) : CondFace(nci) {}
};

// Reflection:
// zero value for normal component and zero gradient for tangential
class CondFaceReflect : public CondFace {
 public:
  CondFaceReflect(size_t nci) : CondFace(nci) {}
};

// Condition for value
template <class V>
class CondFaceVal : public CondFace {
 public:
  CondFaceVal(size_t nci) : CondFace(nci) {}
  virtual V second() const = 0;
};

// Extract single component from a Vect condition
template <class Vect>
class CondFaceValComp : public CondFaceVal<typename Vect::value_type> {
 public:
  using Scal = typename Vect::value_type;
  using P = CondFaceVal<Scal>; // parent
  CondFaceValComp(CondFaceVal<Vect>* o, size_t d)
      : P(o->GetNci()), o_(o), d_(d) {}
  Scal second() const override {
    return o_->second()[d_];
  }

 private:
  CondFaceVal<Vect>* o_;
  size_t d_;
};

// Given value
template <class V>
class CondFaceValFixed : public CondFaceVal<V> {
 public:
  CondFaceValFixed(const V& v, size_t nci) : CondFaceVal<V>(nci), v_(v) {}
  V second() const override {
    return v_;
  }
  void Set(const V& v) {
    v_ = v;
  }

 private:
  V v_;
};

// Condition for gradient
template <class V>
class CondFaceGrad : public CondFace {
 public:
  CondFaceGrad(size_t nci) : CondFace(nci) {}
  virtual V GetGrad() const = 0;
};

// Extract single component from a Vect condition
template <class Vect>
class CondFaceGradComp : public CondFaceGrad<typename Vect::value_type> {
 public:
  using Scal = typename Vect::value_type;
  using P = CondFaceGrad<Scal>; // parent
  CondFaceGradComp(CondFaceGrad<Vect>* o, size_t d)
      : P(o->GetNci()), o_(o), d_(d) {}
  Scal GetGrad() const override {
    return o_->GetGrad()[d_];
  }

 private:
  CondFaceGrad<Vect>* o_;
  size_t d_;
};

// Given gradient
template <class V>
class CondFaceGradFixed : public CondFaceGrad<V> {
 public:
  explicit CondFaceGradFixed(const V& v, size_t nci)
      : CondFaceGrad<V>(nci), v_(v) {}
  V GetGrad() const override {
    return v_;
  }
  void Set(const V& v) {
    v_ = v;
  }

 private:
  V v_;
};

// Evaluates face condition of given type.
template <class T>
UniquePtr<CondFace> Eval(const UniquePtr<CondFace>& b) {
  if (auto d = b.Get<CondFaceVal<T>>()) {
    return UniquePtr<CondFaceValFixed<T>>(d->second(), d->GetNci());
  }
  if (auto d = b.Get<CondFaceGrad<T>>()) {
    return UniquePtr<CondFaceGradFixed<T>>(d->GetGrad(), d->GetNci());
  }
  if (auto d = b.Get<CondFaceExtrap>()) {
    return UniquePtr<CondFaceExtrap>(d->GetNci());
  }
  if (auto d = b.Get<CondFaceReflect>()) {
    return UniquePtr<CondFaceReflect>(d->GetNci());
  }
  throw std::runtime_error(std::string(__func__) + ": Unknown face condition");
}

// Evaluates component i of vector face condition.
template <class Vect>
UniquePtr<CondFace> EvalComp(const UniquePtr<CondFace>& b, size_t i) {
  using Scal = typename Vect::value_type;
  if (auto d = b.Get<CondFaceVal<Vect>>()) {
    return UniquePtr<CondFaceValFixed<Scal>>(d->second()[i], d->GetNci());
  }
  if (auto d = b.Get<CondFaceGrad<Vect>>()) {
    return UniquePtr<CondFaceGradFixed<Scal>>(d->GetGrad()[i], d->GetNci());
  }
  if (auto d = b.Get<CondFaceExtrap>()) {
    return UniquePtr<CondFaceExtrap>(d->GetNci());
  }
  if (auto d = b.Get<CondFaceReflect>()) {
    return UniquePtr<CondFaceReflect>(d->GetNci());
  }
  throw std::runtime_error(std::string(__func__) + ": Unknown face condition");
}

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

using MapCondFace = MapFace<UniquePtr<CondFace>>;

template <class T, class Map>
MapCondFace GetCondZeroGrad(const Map& mf) {
  MapCondFace r;
  for (auto& it : mf) {
    const IdxFace f = it.first;
    r[f].Set<CondFaceGradFixed<T>>(T(0), it.second->GetNci());
  }
  return r;
}

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
  BCond(BCondType type, size_t nci, T val) : type(type), nci(nci), val(val) {}
  BCond(BCondType type, size_t nci) : type(type), nci(nci) {}
  BCondType type = BCondType::dirichlet;
  size_t nci = 0; // neighbor cell id on faces and 0 on embedded boundaries
  T val = T(0); // field value (dirichlet) or normal gradient (neumann)
};
