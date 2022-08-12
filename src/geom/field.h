// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <algorithm>
#include <cassert>
#include <iosfwd>
#include <memory>
#include <stdexcept>
#include <string>

#include "idx.h"
#include "range.h"
#include "util/logger.h"

// Partial implementation of std::vector<T> without specialization for T=bool.
template <class T>
class Vector {
 public:
  Vector() = default;
  ~Vector() {
    if (owning_) {
      delete[] data_;
    }
  }
  explicit Vector(size_t size)
      : size_(size), data_(new T[size]), owning_(true) {}
  Vector(T* data, size_t size) : size_(size), data_(data), owning_(false) {}
  Vector(size_t size, const T& value)
      : size_(size), data_(new T[size]), owning_(true) {
    std::fill(data_, data_ + size_, value);
  }
  Vector(const Vector& other)
      : size_(other.size_), data_(new T[other.size_]), owning_(true) {
    std::copy(other.data_, other.data_ + size_, data_);
  }
  Vector(Vector&& other)
      : size_(other.size_), data_(other.data_), owning_(other.owning_) {
    other.size_ = 0;
    other.data_ = nullptr;
    other.owning_ = true;
  }
  Vector& operator=(const Vector& other) {
    Vector(other).swap(*this);
    return *this;
  }
  Vector& operator=(Vector&& other) {
    if (owning_) {
      delete[] data_;
    }
    size_ = other.size_;
    data_ = other.data_;
    owning_ = other.owning_;
    other.size_ = 0;
    other.data_ = nullptr;
    other.owning_ = true;
    return *this;
  }
  void swap(Vector& other) noexcept {
    std::swap(other, *this);
  }
  size_t size() const noexcept {
    return size_;
  }
  void resize(size_t newsize) {
    if (newsize != size_) {
      Vector(newsize).swap(*this);
    }
  }
  const T* data() const noexcept {
    return data_;
  }
  T* data() noexcept {
    return data_;
  }
  const T& operator[](size_t i) const noexcept {
    return data_[i];
  }
  T& operator[](size_t i) noexcept {
    return data_[i];
  }
  bool owning() const noexcept {
    return owning_;
  }

 private:
  size_t size_ = 0;
  T* data_ = nullptr;
  bool owning_ = true; // true if memory is managed by the object
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
  GField(Value* data, const Range& range)
      : range_(range), data_(data, size(range)) {}
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
  bool owning() const {
    return data_.owning();
  }
  void Reinit(Value* data, const Range& range) {
    GField(data, range).swap(*this);
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
    fassert(
        halo <= halo_, //
        "required " + std::to_string(halo) + " halo cells for field '" +
            GetName() + "', but only '" + std::to_string(halo_) +
            "' are valid");
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

// Writes one component of vector field to given scalar field.
template <class Vect, class Idx>
void GetComponent(
    const GField<Vect, Idx>& fvect, size_t comp,
    GField<typename Vect::value_type, Idx>& fscal) {
  fassert(comp < Vect::dim);
  fscal.Reinit(fvect.GetRange());
  for (auto i : fvect.GetRange()) {
    fscal[i] = fvect[i][comp];
  }
  fscal.SetHalo(fvect.GetHalo());
}

// Returns component of vector field.
template <class Vect, class Idx>
GField<typename Vect::value_type, Idx> GetComponent(
    const GField<Vect, Idx>& fvect, size_t comp) {
  using Scal = typename Vect::value_type;
  GField<Scal, Idx> fscal(fvect.GetRange());
  GetComponent(fvect, comp, fscal);
  return fscal;
}

// Sets component of vector field.
template <class Vect, class Idx>
void SetComponent(
    GField<Vect, Idx>& fvect, size_t comp,
    const GField<typename Vect::value_type, Idx>& fscal) {
  fassert_equal(fscal.size(), fvect.size());
  for (auto i : fvect.GetRange()) {
    fvect[i][comp] = fscal[i];
  }
  fvect.SetHalo(std::min(fvect.GetHalo(), fscal.GetHalo()));
}

// Sets component of vector field.
template <class Vect, class Idx>
void SetComponent(
    GField<Vect, Idx>& fvect, size_t comp, const typename Vect::value_type& s) {
  for (auto i : fvect.GetRange()) {
    fvect[i][comp] = s;
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
  for (auto i : m.template GetRangeIn<Idx>()) {
    out << u[i] << " ";
  }
}
