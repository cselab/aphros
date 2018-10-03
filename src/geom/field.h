#pragma once

#include <ostream>

#include "idx.h"
#include "range.h"

template <class Value_, class Idx_>
class GField {
 public:
  using Idx = Idx_;
  using Value = Value_;
  using Range = GRange<Idx>;

  // Constructs field for empty range
  GField() : d_(nullptr) {}
  // Constructs field for range r
  explicit GField(const Range& r) : r_(r), d_(new Value[r_.size()]) { }
  // Copy constructor
  GField(const GField& o) : GField(o.r_) {
    for (auto i : r_) {
      (*this)[i] = o[i];
    }
  }
  // Move constructor
  GField(GField&& o) : r_(o.r_), d_(o.d_) {
    o.r_.clear();
    o.d_ = nullptr;
  } 
  // Constructs field for range r and initializes with v
  GField(const Range& r, const Value& v) : GField(r) {
    for (auto i : r_) {
      (*this)[i] = v;
    }
  }
  ~GField() {
    Free();
  }
  // Assignment
  GField& operator=(const GField& o) {
    Reinit(o.r_);
    for (auto i : r_) {
      (*this)[i] = o[i];
    }
    return *this;
  }
  // Move assignment
  GField& operator=(GField&& o) {
    Free();
    r_ = o.r_;
    d_ = o.d_;
    o.r_.clear();
    o.d_ = nullptr;
    return *this;
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
    r_.clear();
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
    out << size_t(i) << " " << u[i] << "\n";
  }
  return out;
}

template <class T, class M, class Idx>
void PrintIn(const GField<T, Idx>& u, const M& m, std::ostream& out) {
  for (auto i : m.template GetIn<Idx>()) {
    out << u[i] << " ";
  }
}

