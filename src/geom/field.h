// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <algorithm>
#include <cassert>
#include <iosfwd>
#include <stdexcept>
#include <string>
#include <memory>

#include "idx.h"
#include "range.h"
#include "util/logger.h"

// Partial implementation of std::vector<T> without specialization for T=bool.
template <class T>
class Vector {
 public:
  Vector() = default;
  ~Vector() = default;
  explicit Vector(size_t size) : size_(size), data_(new T[size]) {}
  Vector(size_t size, const T& value) : size_(size), data_(new T[size]) {
    std::fill(data_.get(), data_.get() + size_, value);
  }
  Vector(const Vector& other) : size_(other.size_), data_(new T[other.size_]) {
    std::copy(other.data_.get(), other.data_.get() + size_, data_.get());
  }
  Vector(Vector&& other) = default;
  Vector& operator=(const Vector& other) {
    auto tmp = Vector(other);
    swap(tmp);
    return *this;
  }
  Vector& operator=(Vector&&) = default;
  void swap(Vector& other) {
    std::swap(other, *this);
  }
  size_t size() const {
    return size_;
  }
  void resize(size_t newsize) {
    if (newsize != size_) {
      auto tmp = Vector(newsize);
      swap(tmp);
    }
  }
  const T* data() const {
    return data_.get();
  }
  T* data() {
    return data_.get();
  }
  const T& operator[](size_t i) const {
    return data_.get()[i];
  }
  T& operator[](size_t i) {
    return data_.get()[i];
  }

 private:
  size_t size_ = 0;
  std::unique_ptr<T[]> data_;
};


template <class Value_, class Idx_>
class GField {
 public:
  using Idx = Idx_;
  using Value = Value_;
  using Range = GRange<Idx>;
  static constexpr int kMaxHalo = 128;

  GField() = default;
  explicit GField(const Range& range) : range_(range), data_(size(range)) {}
  ~GField() = default;
  GField(const GField& o) = default;
  GField(GField&& o) = default;
  GField(const Range& range, const Value& value)
      : range_(range), data_(size(), value) {}
  GField& operator=(const GField& o) = default;
  GField& operator=(GField&& o) = default;
  static size_t size(const Range& range) {
    return static_cast<size_t>(*range.end());
  }
  size_t size() const {
    return size(range_);
  }
  bool empty() const {
    return size() == 0;
  }
  // Changes the range, reallocates memory if size differs
  void Reinit(const Range& range) {
    if (range != range_) {
      range_ = range;
      data_.resize(size(range));
    }
  }
  // Changes the range, reallocates memory if size differs.
  // Overwrites data with `value`.
  void Reinit(const Range& range, const Value& value) {
    Reinit(range);
    for (auto i : range_) {
      (*this)[i] = value;
    }
  }
  Range GetRange() const {
    return range_;
  }
  void swap(GField& o) {
    std::swap(o, *this);
  }
  Value* data() {
    return data_.data();
  }
  const Value* data() const {
    return data_.data();
  }
  Value& operator[](Idx i) {
    assert(size_t(i) < size());
    return data_[size_t(i)];
  }
  const Value& operator[](Idx i) const {
    assert(size_t(i) < size());
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
  Vector<Value> data_;
  std::string name_;
  int halo_ = kMaxHalo;
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
