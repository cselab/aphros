#pragma once

#include <ostream>

#include "range.h"

template <class Value_, class Idx_>
class GField {
 public:
  using Idx = Idx_;
  using Value = Value_;
  using Range = GRange<Idx>;

  // Constructs field for range of zero size
  GField() : d_(nullptr) {}

  // Constructs field for range r
  explicit GField(const Range& r) : r_(r), d_(new Value[r_.size()]) {}

  // Constructs field for range r and initializes with v
  GField(const Range& r, const Value& v) : GField(r) {
    for (auto i : r_) {
      (*this)[i] = v;
    }
  }
  ~GField() {
    Free();
  }
  // Copies field of different type
  template <class U>
  GField(const GField<U, Idx>& o) : GField(o.r_) {
    for (auto i : r_) {
      (*this)[i] = o[i];
    }
  }
  size_t size() const { 
    return r_.size(); 
  }
  bool empty() const {
    return size() == 0;
  }
  // Changes range to r, reallocates memory if size differs
  void Reinit(const Range& r) {
    if (r != r_) {
      GField(r).swap(*this);
    }
  }
  // Changes range to r, reallocates memory if size differs, overwrites with v
  void Reinit(const Range& r, const Value& v) {
    Reinit(r);
    for (auto i : r_) {
      (*this)[i] = v;
    }
  }
  // Deallocates memory, resets range to zero size
  void Free() {
    Range().swap(r_);
    delete[] d_;
    d_ = nullptr;
  }
  Range GetRange() const {
    return r_;
  }
  void swap(GField& o) {
    r_.swap(o.r_);
    std::swap(d_, o.d_);
  }
  Value* data() {
    return d_;
  }
  const Value* data() const {
    return d_;
  }
  Value& operator[](const Idx& i) {
    assert(size_t(i) >= 0 && size_t(i) < r_.size());
    return d_[size_t(i)];
  }
  const Value& operator[](const Idx& i) const {
    assert(size_t(i) >= 0 && size_t(i) < r_.size());
    return d_[size_t(i)];
  }

 private:
  template <class, class>
  friend class GField;

  Range r_;
  Value* d_;
};

template <class T>
using FieldCell = GField<T, IdxCell>;

template <class T>
using FieldFace = GField<T, IdxFace>;

template <class T>
using FieldNode = GField<T, IdxNode>;

// Write component of vector field 
// to given scalar field (resize if needed).
// Requires: defined Vect::value_type.
template <class Vect, class Idx>
void GetComponent(
    const GField<Vect, Idx>& fv, size_t n, 
    GField<typename Vect::value_type, Idx>& fs) {
  fs.Reinit(fv.GetRange());
  for (auto i : fv.GetRange()) {
    fs[i] = fv[i][n];
  }
}

// Return component of vector field.
// Requires: defined Vect::value_type.
template <class Vect, class Idx>
GField<typename Vect::value_type, Idx> GetComponent(
    const GField<Vect, Idx>& fv, size_t n) {
  using Scal = typename Vect::value_type;
  GField<Scal, Idx> fs(fv.GetRange());
  GetComponent(fv, n, fs);
  return fs;
}

// Set component of vector field.
// Requires: defined Vect::value_type.
template <class Vect, class Idx>
void SetComponent(
    GField<Vect, Idx>& fv, size_t n, 
    const GField<typename Vect::value_type, Idx>& fs) {
  for (auto i : fv.GetRange()) {
    fv[i][n] = fs[i];
  }
}

template <class T, class Idx>
std::ostream& operator<<(std::ostream& out, const GField<T, Idx>& u) {
  for (auto i : u.GetRange()) {
    out << i.GetRaw() << " " << u[i] << "\n";
  }
  return out;
}



