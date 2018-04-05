#pragma once

#include <vector>
#include <ostream>

#include "range.h"

namespace geom {

template <class Value_, class Idx_>
class GField {
 public:
  using Idx = Idx_;
  using Value = Value_;
  using Range = GRange<Idx>;
  using Cont = std::vector<Value>; // container

  GField() {}
  template <class U>
  GField(const GField<U, Idx>& o /*other*/)
    : d_(o.d_.begin(), o.d_.end()) {}
  explicit GField(const Range& r)
    : d_(r.size()) {}
  GField(const Range& r, const Value& val)
    : d_(r.size(), val) {}
  size_t size() const {
    return d_.size();
  }
  void Reinit(const Range& r) {
    d_.resize(r.size());
  }
  void Reinit(const Range& r, const Value& v) {
    d_.assign(r.size(), v);
  }
  Range GetRange() const {
    return Range(0, size());
  }
  void resize(size_t size) {
    d_.resize(size);
  }
  bool empty() const {
    return d_.empty();
  }
  // Cont::pointer (instead of Value*) required for vector<bool>
  typename Cont::pointer data() {
    return d_.data();
  }
  typename Cont::const_pointer data() const {
    return d_.data();
  }
  void push_back(const Value& v) {
    d_.push_back(v);
  }
  typename Cont::reference operator[](const Idx& idx) {
    return d_[idx.GetRaw()];
  }
  typename Cont::const_reference operator[](const Idx& idx) const {
    return d_[idx.GetRaw()];
  }

 private:
  // Required for constructor copying from another type
  template <class OValue, class OIdx>
  friend class GField;

  std::vector<Value> d_;
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


} // namespace geom
