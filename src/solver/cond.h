#pragma once

#include <memory>
#include <utility>

#include "geom/map.h"

namespace solver {

class CondFace {
 public:
  virtual ~CondFace() {}
  virtual size_t GetNci() const { return nci_; }

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
  virtual V GetValue() const = 0;
};

// Extract single component from a Vect condition
template <class Vect>
class CondFaceValComp : public CondFaceVal<typename Vect::value_type> {
 public:
  using Scal = typename Vect::value_type;
  using P = CondFaceVal<Scal>; // parent
  CondFaceValComp(CondFaceVal<Vect>* o, size_t d)
      : P(o->GetNci()), o_(o), d_(d) {}
  Scal GetValue() const override { return o_->GetValue()[d_]; }

 private:
  CondFaceVal<Vect>* o_;
  size_t d_;
};

// Given value
template <class V>
class CondFaceValFixed : public CondFaceVal<V> {
 public:
  CondFaceValFixed(const V& v, size_t nci) : CondFaceVal<V>(nci), v_(v) {}
  V GetValue() const override { return v_; }
  void Set(const V& v) { v_ = v; }

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
  Scal GetGrad() const override { return o_->GetGrad()[d_]; }

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
  V GetGrad() const override { return v_; }
  void Set(const V& v) { v_ = v; }

 private:
  V v_;
};

class CondCell {
 public:
  virtual ~CondCell() {}
};

// Condition for value
template <class V>
class CondCellVal : public CondCell {
 public:
  virtual V GetValue() const = 0;
};

// Given value
template <class V>
class CondCellValFixed : public CondCellVal<V> {
 public:
  explicit CondCellValFixed(const V& v) : v_(v) {}
  V GetValue() const override { return v_; }
  void Set(const V& v) { v_ = v; }

 private:
  V v_;
};

} // namespace solver

template <class T>
class UniquePtr {
 public:
  UniquePtr() = default;
  UniquePtr(const UniquePtr&) = delete;
  UniquePtr(UniquePtr&&) = default;
  template <class U>
  UniquePtr(const UniquePtr<U>&& o) 
      : p_(o.p_) {}
  UniquePtr(std::nullptr_t) {}
  UniquePtr& operator=(const UniquePtr&) = delete;
  UniquePtr& operator=(UniquePtr&&) = default;
  template <class U, class ... Args>
  UniquePtr(Args ... args) 
      : p_(new  U(std::forward<Args>(args)...)) {}
  template <class U, class ... Args>
  void Set(Args ... args) {
    p_ = std::unique_ptr<T>(new U(std::forward<Args>(args)...));
  }
  void Set(std::nullptr_t) {
    p_.reset(nullptr);
  }
  void Set(int) = delete;
  template <class U=T>
  U* Get() {
    return dynamic_cast<U*>(p_.get());
  }
  template <class U=T>
  const U* Get() const {
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
};

using MapCondFace = MapFace<UniquePtr<solver::CondFace>>;
using MapCondFaceFluid = MapFace<UniquePtr<solver::CondFace>>;
