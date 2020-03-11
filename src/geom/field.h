// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <cassert>
#include <iosfwd>
#include <stdexcept>
#include <string>

#include "idx.h"
#include "range.h"
#include "util/logger.h"

template <class Value_, class Idx_>
class GField {
 public:
  using Idx = Idx_;
  using Value = Value_;
  using Range = GRange<Idx>;

  // Constructs field for empty range
  GField() : data_(nullptr) {}
  // Constructs field for range r
  explicit GField(const Range& r) : range_(r), data_(new Value[size()]) {}
  // Copy constructor
  GField(const GField& o) : GField(o.range_) {
    for (auto i : range_) {
      (*this)[i] = o[i];
    }
    name_ = o.name_;
    halo_ = o.halo_;
  }
  // Move constructor
  GField(GField&& o)
      : range_(o.range_), data_(o.data_), name_(o.name_), halo_(o.halo_) {
    o.range_.clear();
    o.data_ = nullptr;
    o.name_.clear();
  }
  // Constructs field for range r and initializes with v
  GField(const Range& r, const Value& v) : GField(r) {
    for (auto i : range_) {
      (*this)[i] = v;
    }
  }
  ~GField() {
    Free();
  }
  // Assignment
  GField& operator=(const GField& o) {
    if (this != &o) {
      Reinit(o.range_);
      for (auto i : range_) {
        (*this)[i] = o[i];
      }
      name_ = o.name_;
      halo_ = o.halo_;
    }
    return *this;
  }
  // Move assignment
  GField& operator=(GField&& o) {
    Free();
    range_ = o.range_;
    data_ = o.data_;
    name_ = o.name_;
    halo_ = o.halo_;
    o.range_.clear();
    o.data_ = nullptr;
    o.name_.clear();
    return *this;
  }
  size_t size() const {
    return static_cast<size_t>(*range_.end());
  }
  bool empty() const {
    return size() == 0;
  }
  // Changes range to r, reallocates memory if size differs
  void Reinit(const Range& r) {
    if (r != range_) {
      GField(r).swap(*this);
    }
  }
  // Changes range to r, reallocates memory if size differs, overwrites with v
  void Reinit(const Range& r, const Value& v) {
    Reinit(r);
    for (auto i : range_) {
      (*this)[i] = v;
    }
  }
  // Deallocates memory, resets range to zero size
  void Free() {
    range_.clear();
    delete[] data_;
    data_ = nullptr;
  }
  Range GetRange() const {
    return range_;
  }
  void swap(GField& o) {
    std::swap(range_, o.range_);
    std::swap(data_, o.data_);
    std::swap(name_, o.name_);
    std::swap(halo_, o.halo_);
  }
  Value* data() {
    return data_;
  }
  const Value* data() const {
    return data_;
  }
  Value& operator[](const Idx& i) {
    assert(size_t(i) >= 0 && size_t(i) < size_t(*range_.end()));
    return data_[size_t(i)];
  }
  const Value& operator[](const Idx& i) const {
    assert(size_t(i) >= 0 && size_t(i) < size_t(*range_.end()));
    return data_[size_t(i)];
  }
  std::string GetName() const {
    return name_;
  }
  void SetName(const std::string name) {
    name_ = name;
  }
  void CheckHalo(int halo) const {
    if (halo_ < halo) {
      throw std::runtime_error(
          FILELINE + ": required " + std::to_string(halo) +
          " halos for field '" + GetName() + "' but only " +
          std::to_string(halo_) + " are valid");
    }
  }
  void SetHalo(int halo) {
    halo_ = halo;
  }
  void LimitHalo(int halo) {
    if (halo < halo_) {
      halo_ = halo;
    }
  }
  int GetHalo() const {
    return halo_;
  }

 private:
  Range range_;
  Value* data_;
  std::string name_;
  int halo_ = 100; // default to max
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
  fs.SetHalo(fv.GetHalo());
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
  fv.SetHalo(std::min(fv.GetHalo(), fs.GetHalo()));
}

// Set component of vector field.
// Requires: defined Vect::value_type.
template <class Vect, class Idx>
void SetComponent(
    GField<Vect, Idx>& fv, size_t n, const typename Vect::value_type& s) {
  for (auto i : fv.GetRange()) {
    fv[i][n] = s;
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
